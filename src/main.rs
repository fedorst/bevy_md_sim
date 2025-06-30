mod config;
mod core;
mod debug;
mod setup;
mod simulation;
mod ui;
mod visualization;

use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCameraPlugin;
use clap::Parser;
use core::CorePlugin;
use core::MoleculeSelection; // Import our new resource
use core::SimulationParameters;
use core::SimulationState;
use core::Thermostat;
use debug::DebugPlugin;
use setup::SetupPlugin;
use simulation::SimulationPlugin;
use ui::UIPlugin;
use visualization::VisualizationPlugin;

// | ID | Feature | Status | Effort | Notes |
// | :--- | :--- | :--- | :--- | :--- |
// | 14 | **Display Total System Energy** | 💡 **NEXT UP** | **LOW** | Huge value for simulation stability! A new system to sum kinetic & potential energy is straightforward and can be displayed in a new UI text element. |
// | 3 | **Delete Selected Atom** | 💡 **NEXT UP** | **MEDIUM** | Very satisfying interactive feature. The core is just `commands.despawn()`. The "medium" effort comes from needing to carefully remove the atom's bonds/angles from `SystemConnectivity` to prevent crashes. |
// | --- | | | | |
// | 7 | **Create/Delete Bonds (Shift+Select)**| To Do | **MEDIUM** | Involves managing multi-selection state, checking for key modifiers, and editing `SystemConnectivity`. A great next step for interactivity. |
// | 9 | **Implement Double Bonds** | To Do | **MEDIUM** | This is more of a core simulation feature. It requires changing `Bond` data structures and the force calculation logic. |
// | 13 | **Further UI Improvements** | To Do | **MEDIUM** | A good example would be adding buttons for "Pause/Play" or "Reset Simulation" instead of using keyboard-only controls. |
// | --- | | | | |
// | 1 | **Marquee (Rectangle) Select** | To Do | **HIGH** | Requires complex logic for 2D-to-3D space conversion, managing a selection box, and checking many entities against it. |
// | 4 | **Draggable Atoms** | To Do | **HIGH** | Requires implementing a 3D translation gizmo or integrating a library for it, which is a significant task. |
// | 6 | **Add Atoms from a Sidebar** | To Do | **HIGH** | This is a major feature requiring significant UI work (e.g., `bevy_egui`) and new spawning/state management logic. |
// | 11| **Playback System** | To Do | **HIGH** | Involves serializing the entire system state each frame and building a separate playback mode and UI. Very complex. |

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
            DebugPlugin,
        ))
        .run();
}
