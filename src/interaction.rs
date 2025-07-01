use crate::components::{AtomType, BondVisualization, Force, Velocity};
use crate::resources::SystemConnectivity;
use crate::simulation::PhysicsSet;
use bevy::color::palettes::basic::{BLACK, BLUE, RED};
use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCamera;
use bevy_picking::prelude::MeshPickingPlugin;
use bevy_picking::{
    events::{Click, Pointer},
    hover::PickingInteraction,
};

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct InteractionSet;

pub struct InteractionPlugin;

impl Plugin for InteractionPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(MeshPickingPlugin)
            .init_resource::<SelectionState>()
            .add_observer(handle_selection)
            .add_systems(
                Update,
                (
                    update_highlights,
                    draw_selection_gizmos,
                    delete_selected_atom,
                    focus_camera_on_selection,
                )
                    .in_set(InteractionSet)
                    .after(PhysicsSet),
            );
    }
}

#[derive(Resource, Default)]
pub struct SelectionState {
    pub selected: Vec<Entity>,
}

fn handle_selection(
    trigger: Trigger<Pointer<Click>>,
    mut selection: ResMut<SelectionState>,
    atoms: Query<&AtomType>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    let target = trigger.target();
    if atoms.get(target).is_ok() {
        let shift_pressed = keys.pressed(KeyCode::ShiftLeft) || keys.pressed(KeyCode::ShiftRight);
        if shift_pressed {
            if let Some(index) = selection.selected.iter().position(|&e| e == target) {
                selection.selected.remove(index);
                info!("Deselected {:?} from multi-select.", target);
            } else {
                selection.selected.push(target);
                info!("Added {:?} to multi-select.", target);
            }
        } else {
            if selection.selected.len() == 1 && selection.selected[0] == target {
                selection.selected.clear();
                info!("Deselected all.");
            } else {
                selection.selected.clear();
                selection.selected.push(target);
                info!("Selected only {:?}.", target);
            }
        }
    }
}

fn delete_selected_atom(
    mut commands: Commands,
    keys: Res<ButtonInput<KeyCode>>,
    mut selection: ResMut<SelectionState>,
    mut connectivity: ResMut<SystemConnectivity>,
    bond_vis_query: Query<(Entity, &BondVisualization)>,
) {
    if keys.just_pressed(KeyCode::Delete) {
        let entities_to_delete = selection.selected.clone();
        if !entities_to_delete.is_empty() {
            info!(
                "Deletion key pressed for entities: {:?}",
                entities_to_delete
            );

            for entity_to_delete in &entities_to_delete {
                commands.entity(*entity_to_delete).despawn();

                connectivity.bonds.retain(|bond| {
                    if &bond.a == entity_to_delete || &bond.b == entity_to_delete {
                        for (vis_entity, bond_vis) in &bond_vis_query {
                            if (bond_vis.atom1 == bond.a && bond_vis.atom2 == bond.b)
                                || (bond_vis.atom1 == bond.b && bond_vis.atom2 == bond.a)
                            {
                                commands.entity(vis_entity).despawn();
                                info!(
                                    "Despawned bond visualization for bond between {:?} and {:?}",
                                    bond.a, bond.b
                                );
                            }
                        }
                        return false;
                    }
                    // Keep the bond if it's not connected to the deleted atom.
                    true
                });
                connectivity.angles.retain(|angle| {
                    &angle.a != entity_to_delete
                        && &angle.b != entity_to_delete
                        && &angle.center != entity_to_delete
                });
            }

            selection.selected.clear(); // Clear the entire selection
            info!("Cleared selection.");
        }
    }
}

fn update_highlights(
    selection: Res<SelectionState>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    query: Query<
        (
            Entity,
            &PickingInteraction,
            &MeshMaterial3d<StandardMaterial>,
        ),
        With<AtomType>,
    >,
) {
    for (entity, interaction, material_handle_comp) in &query {
        if let Some(material) = materials.get_mut(&material_handle_comp.0) {
            material.emissive = BLACK.into();

            if selection.selected.contains(&entity) {
                material.emissive = material.emissive.lighter(0.3);
            } else {
                if *interaction == PickingInteraction::Hovered {
                    material.emissive = material.emissive.lighter(0.1);
                }
            }
        }
    }
}

fn draw_selection_gizmos(
    selection: Res<SelectionState>,
    atom_query: Query<(&Transform, &Force, &Velocity), With<AtomType>>,
    mut gizmos: Gizmos,
) {
    gizmos.clear();

    for &selected_entity in &selection.selected {
        // The rest of the logic is the same, but it now runs for each selected atom.
        if let Ok((transform, force, velocity)) = atom_query.get(selected_entity) {
            let start = transform.translation;
            let force_vector = force.total * 0.0004;
            let f_end = start + force_vector;
            let v_end = start + velocity.0 * 0.1;
            gizmos.arrow(start, f_end, RED);
            gizmos.arrow(start, v_end, BLUE);
        }
    }
}

fn focus_camera_on_selection(
    keys: Res<ButtonInput<KeyCode>>,
    selection: Res<SelectionState>,
    mut camera_query: Query<&mut PanOrbitCamera>,
    atom_query: Query<&Transform, With<AtomType>>,
) {
    if keys.just_pressed(KeyCode::KeyF) {
        if let Some(&entity_to_focus_on) = selection.selected.last() {
            if let Ok(mut pan_orbit_camera) = camera_query.single_mut() {
                // The query's .get() method takes an owned Entity, which `entity_to_focus_on` now is.
                if let Ok(atom_transform) = atom_query.get(entity_to_focus_on) {
                    info!("Focusing camera on {:?}", entity_to_focus_on);
                    pan_orbit_camera.target_focus = atom_transform.translation;
                }
            }
        }
    }
}
