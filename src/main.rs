mod core;
mod debug;
mod setup;
mod simulation;
mod ui;
mod visualization;

use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCameraPlugin;
use core::CorePlugin;
use debug::DebugPlugin;
use setup::SetupPlugin;
use simulation::SimulationPlugin;
use ui::UIPlugin;
use visualization::VisualizationPlugin;

// | ID | Feature | Status | Effort | Notes |
// | :--- | :--- | :--- | :--- | :--- |
// | 10 | **Add Nitrogen Atom** | 💡 **NEXT UP** | **LOW** | Easy to implement. Just requires adding to the `AtomType` enum and finding/adding its force field parameters. Great "quick win". |
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

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
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

/*
fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .init_resource::<ForceField>() // Add our constants as a resource
        .init_resource::<SystemConnectivity>()
        .init_resource::<StepSimulation>()
        .insert_resource(SimulationParameters { dt: 1e-4 })
        .init_resource::<StepCount>()
        .init_resource::<SimulationState>()
        .add_systems(Startup, (setup, build_exclusions, run_once).chain())
        .add_systems(
            Update,
            (
                // debug_state,
                // --- VERLET STEP 1: Update positions and half-step velocity ---
                integrate_position_and_half_velocity,
                // --- VERLET STEP 2: Calculate new forces based on new positions ---
                reset_forces,
                calculate_bond_forces,
                calculate_angle_forces,
                calculate_non_bonded_forces,
                // log_total_forces,
                // --- VERLET STEP 3: Complete the velocity update with new forces ---
                finish_velocity_update,
                // integrate_position,
                // update_velocity,
                update_bond_visuals,
                end_simulation_step,
            )
                .chain()
                .run_if(resource_equals(StepSimulation(true))),
        )
        .add_systems(
            Update,
            (
                handle_simulation_control,
                continuous_simulation,
                update_pause_text,
            ),
        )
        .run();
}

fn debug_state(
    query: Query<(Entity, &AtomType, &Transform)>,
    connectivity: Res<SystemConnectivity>,
) {
    info!("--- CURRENT STATE ---");
    // This is a bit more complex now, but it's the correct "Bevy way".
    // We need to find the entities from our connectivity resource.
    let o_entity = connectivity.bonds[0].a;
    let h1_entity = connectivity.bonds[0].b;
    let h2_entity = connectivity.bonds[1].b;

    // Use the query to get the transform data.
    let o_pos = query.get(o_entity).unwrap().2.translation;
    let h1_pos = query.get(h1_entity).unwrap().2.translation;
    let h2_pos = query.get(h2_entity).unwrap().2.translation;

    info!("Oxygen at: {:?}", o_pos);
    info!("Hydrogen 1 at: {:?}", h1_pos);
    info!("Hydrogen 2 at: {:?}", h2_pos);
    info!("Distance O-H1: {}", o_pos.distance(h1_pos));
    info!("Distance O-H2: {}", o_pos.distance(h2_pos));
    info!("---------------------------------\n");
}

fn log_total_forces(query: Query<(&Force, &AtomType)>) {
    let mut total_force = Vec3::ZERO;
    for (force, atom_type) in &query {
        total_force += force.0;
        info!(
            "Atom {:?}: F={:?} (|F|={:.6})",
            atom_type,
            force.0,
            force.0.length()
        );
    }
    info!("TOTAL SYSTEM FORCE: {:?}", total_force);
}
 */
