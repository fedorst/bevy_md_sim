mod components;
mod config;
mod cursor;
mod interaction;
mod persistence;
mod resources;
mod setup;
mod simulation;
mod spawning;
mod ui;
mod visualization;

use bevy::prelude::*;
use bevy_egui::EguiPlugin;
use bevy_panorbit_camera::PanOrbitCameraPlugin;
use clap::Parser;
use cursor::CustomCursorPlugin;
use interaction::{InteractionPlugin, RebuildConnectivityEvent};
use persistence::PersistencePlugin;
use resources::*;
use setup::SetupPlugin;
use simulation::SimulationPlugin;
use spawning::SpawningPlugin;
use ui::UIPlugin;
use visualization::VisualizationPlugin;

#[derive(States, Debug, Clone, PartialEq, Eq, Hash, Default)]
pub enum AppState {
    #[default]
    Initializing, // Start here to allow setup systems to run.
    Running,
    Paused,
}

#[derive(Resource, Parser, Debug)]
#[command(version, about, long_about = None)]
struct CliArgs {
    /// The name of the molecule to load from the assets/molecules/{molecule}.json.
    #[arg(short, long, default_value = "ethanol")]
    molecule: String,

    /// step length in picoseconds.
    #[arg(long, default_value_t = 1e-4)]
    dt: f32,

    /// Target temperature in Kelvin.
    #[arg(long, default_value_t = 300.0)]
    temp: f32,

    /// Thermostat coupling time constant (tau) in picoseconds.
    #[arg(long, default_value_t = 0.001)]
    tau: f32,

    /// Start the simulation in a paused state.
    #[arg(long, default_value_t = false)]
    paused: bool,
}

pub struct CorePlugin;

impl Plugin for CorePlugin {
    fn build(&self, app: &mut App) {
        let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
        let ff_path = std::path::Path::new(&manifest_dir)
            .join("assets/force_field.json")
            .to_str()
            .unwrap()
            .to_string();

        app.insert_resource(ForceField::from_file(&ff_path))
            .insert_resource::<ForceMultiplier>(ForceMultiplier(1.0))
            .init_resource::<SystemConnectivity>()
            .init_resource::<StepCount>()
            .init_resource::<PauseMenuState>()
            .init_resource::<ActiveWallTime>()
            .init_resource::<ThermostatScale>()
            .init_resource::<SystemEnergy>()
            .init_resource::<EnergyHistory>()
            .init_resource::<CurrentTemperature>()
            .init_resource::<AtomIdMap>()
            .init_resource::<LastSaveTime>()
            .init_resource::<LastClick>()
            .init_resource::<ExcludedPairs>();
        // .add_systems(Update, track_active_wall_time);
    }
}

fn check_initialization_complete(
    mut commands: Commands,
    cli: Res<CliArgs>,
    // This system runs when a RebuildConnectivityEvent is sent, which is our signal that setup is done.
    mut rebuild_reader: EventReader<RebuildConnectivityEvent>,
    mut next_state: ResMut<NextState<AppState>>,
) {
    if rebuild_reader.read().last().is_some() {
        info!("Initialization complete. Transitioning to target state.");
        if cli.paused {
            next_state.set(AppState::Paused);
        } else {
            next_state.set(AppState::Running);
        }
        // This system has done its job, so we despawn it.
        commands.remove_resource::<CliArgs>();
    }
}

fn track_active_wall_time(
    time: Res<Time>,
    app_state: Res<State<AppState>>,
    mut active_wall_time: ResMut<ActiveWallTime>,
) {
    if *app_state.get() == AppState::Running {
        active_wall_time.0 += time.delta_secs();
    }
}

fn main() {
    let args = CliArgs::parse();
    info!("CLI arguments parsed. Loading molecule: {}", args.molecule);

    App::new()
        .add_plugins(DefaultPlugins)
        .init_state::<AppState>()
        .insert_resource(SimulationParameters { dt: args.dt })
        .insert_resource(Thermostat {
            target_temperature: args.temp,
            tau: args.tau,
        })
        .insert_resource(args)
        .add_plugins((
            PanOrbitCameraPlugin,
            EguiPlugin::default(),
            CorePlugin,
            SetupPlugin,
            SimulationPlugin,
            UIPlugin,
            VisualizationPlugin,
            InteractionPlugin,
            SpawningPlugin,
            CustomCursorPlugin,
            PersistencePlugin,
        ))
        .add_systems(
            Update,
            (
                check_initialization_complete.run_if(in_state(AppState::Initializing)),
                track_active_wall_time,
            ),
        )
        .run();
}
