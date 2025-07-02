mod components;
mod config;
mod interaction;
mod resources;
mod setup;
mod simulation;
mod ui;
mod visualization;

use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCameraPlugin;
use clap::Parser;
use interaction::InteractionPlugin;
use resources::*;
use setup::SetupPlugin;
use simulation::SimulationPlugin;
use ui::UIPlugin;
use visualization::VisualizationPlugin;

#[derive(Parser, Debug)]
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
            .init_resource::<SystemConnectivity>()
            .init_resource::<StepSimulation>()
            .init_resource::<StepCount>()
            .init_resource::<PauseMenuState>()
            .init_resource::<SimulationState>()
            .init_resource::<ActiveWallTime>()
            .init_resource::<ThermostatScale>()
            .init_resource::<SystemEnergy>()
            .init_resource::<CurrentTemperature>()
            .init_resource::<AtomIdMap>()
            .init_resource::<ExcludedPairs>();
    }
}

fn main() {
    let args = CliArgs::parse();
    info!("CLI arguments parsed. Loading molecule: {}", args.molecule);

    App::new()
        .add_plugins(DefaultPlugins)
        .insert_resource(MoleculeSelection(args.molecule))
        .insert_resource(SimulationParameters { dt: args.dt })
        .insert_resource(Thermostat {
            target_temperature: args.temp,
            tau: args.tau,
        })
        .insert_resource(SimulationState {
            paused: args.paused,
        })
        .add_plugins((
            PanOrbitCameraPlugin,
            CorePlugin,
            SetupPlugin,
            SimulationPlugin,
            UIPlugin,
            VisualizationPlugin,
            InteractionPlugin,
        ))
        .run();
}
