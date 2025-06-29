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

/*
 * further ideas
 * 1. select atoms by dragging across the viewport (green rectangle). mabye some examples exist.
 * 2. selected atoms display their forces using https://bevy.org/examples/gizmos/axes/ : aggregate force direction and magnitude. and also same for velocity? in a different color
 * 3. can press del to delete selected atoms
 * 4. if a single atom is selected, it is draggable using axis helpers
 * 5. proper orbital camera is a requirement then, i guess. there should be some examples...
 * 6. in a sidebar, i can add new atoms to the scene. and then drag them
 * 7. if i shift+select 2 atoms, i should be able to create a bond by pressing b. if one already exists between them, pressing b will remove it
 * 8. also make it possible to have a carbon atom
 * 9. also make it possible to have a double bond
 * 10. Also make it possible to have a nitrogen atom
 * 11. make playback of already-simulated stuff possible. with a play button or sth. record all intermediate states per each frame.
 * 12. i understand that not all of this is immediately doable. lets focus on what creates most value with least effort and prioritize our efforts accordingly.
 * 13. further UI improvements
 * 14. display total force or enthalpy of selected system.
 */

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
