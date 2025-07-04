use crate::components::{Atom, BondVisualization, Force, Velocity};
use crate::resources::{
    Angle, Bond, BondOrder, ForceField, SharedAssetHandles, SystemConnectivity,
};
use crate::simulation::PhysicsSet;
use crate::visualization::VisualizationSet;
use bevy::color::palettes::basic::{BLACK, BLUE, RED};
use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCamera;
use bevy_picking::prelude::MeshPickingPlugin;
use bevy_picking::{
    events::{Click, Pointer},
    hover::PickingInteraction,
};
use std::collections::HashMap;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct InteractionSet;

#[derive(Event, Debug)]
pub struct RebuildConnectivityEvent;

pub struct InteractionPlugin;

impl Plugin for InteractionPlugin {
    fn build(&self, app: &mut App) {
        app.configure_sets(Update, VisualizationSet.after(InteractionSet))
            .add_event::<RebuildConnectivityEvent>()
            .add_plugins(MeshPickingPlugin)
            .init_resource::<SelectionState>()
            .add_observer(handle_selection)
            .add_systems(
                Update,
                (
                    (manage_bonds, delete_selected_atom),
                    ApplyDeferred,
                    (
                        rebuild_on_event,
                        update_highlights,
                        draw_selection_gizmos,
                        focus_camera_on_selection,
                    ),
                )
                    .chain()
                    .in_set(InteractionSet)
                    .after(PhysicsSet),
            );
    }
}
fn manage_bonds(
    keys: Res<ButtonInput<KeyCode>>,
    selection: Res<SelectionState>,
    mut connectivity: ResMut<SystemConnectivity>,
    force_field: Res<ForceField>,
    atom_query: Query<&Atom>,
    // It now takes an EventWriter to signal that a rebuild is needed.
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
) {
    if keys.just_pressed(KeyCode::KeyB) && selection.selected.len() == 2 {
        let e1 = selection.selected[0];
        let e2 = selection.selected[1];

        let bond_index = connectivity
            .bonds
            .iter()
            .position(|b| (b.a == e1 && b.b == e2) || (b.a == e2 && b.b == e1));

        if let Some(index) = bond_index {
            info!("Deleting bond between {:?} and {:?}", e1, e2);
            connectivity.bonds.remove(index);
        } else {
            let Ok(type1) = atom_query.get(e1) else {
                return;
            };
            let Ok(type2) = atom_query.get(e2) else {
                return;
            };

            if force_field.bond_params.contains_key(&(
                type1.type_name.clone(),
                type2.type_name.clone(),
                BondOrder::Single,
            )) {
                info!("Creating bond between {:?} and {:?}", e1, e2);
                connectivity.bonds.push(Bond {
                    a: e1,
                    b: e2,
                    order: BondOrder::Single,
                });
            } else {
                warn!(
                    "Cannot create bond: No parameters found for bond type {:?}-{:?}",
                    type1.type_name, type2.type_name
                );
                return;
            }
        }

        // Instead of calling a function, just send an event.
        rebuild_writer.write(RebuildConnectivityEvent);
    }
}

fn rebuild_on_event(
    mut rebuild_reader: EventReader<RebuildConnectivityEvent>,
    mut commands: Commands,
    mut connectivity: ResMut<SystemConnectivity>,
    bond_vis_query: Query<(Entity, &BondVisualization)>,
    shared_handles: Res<SharedAssetHandles>,
    force_field: Res<ForceField>,
    atom_query: Query<&Atom>,
) {
    // Only run if the event was actually sent. This is an efficient way to check.
    if rebuild_reader.is_empty() {
        return;
    }
    // Important: Consume the events to prevent this from running again next frame.
    rebuild_reader.clear();

    info!("Rebuilding angles and bond visuals...");

    // 1. Despawn all existing bond visuals to ensure a clean slate
    for (entity, _) in bond_vis_query.iter() {
        commands.entity(entity).despawn();
    }

    // 2. Clear all existing angles, they will be regenerated
    connectivity.angles.clear();

    // 3. Re-spawn visuals and regenerate angles from the new bond list
    let mut adjacency: HashMap<Entity, Vec<Entity>> = HashMap::new();
    for bond in &connectivity.bonds {
        // Build adjacency list for angle generation
        adjacency.entry(bond.a).or_default().push(bond.b);
        adjacency.entry(bond.b).or_default().push(bond.a);

        // Check if this bond is defined in the force field
        let Ok([atom1, atom2]) = atom_query.get_many([bond.a, bond.b]) else {
            continue;
        };
        let key = (atom1.type_name.clone(), atom2.type_name.clone(), bond.order);

        let material_handle = if force_field.bond_params.contains_key(&key) {
            shared_handles.bond_material.clone()
        } else {
            warn!(
                "Found undefined bond between types: {} and {}",
                atom1.type_name, atom2.type_name
            );
            shared_handles.undefined_bond_material.clone()
        };

        let num_strands = match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        };
        for i in 0..num_strands {
            commands.spawn((
                Mesh3d(shared_handles.bond_mesh.clone()),
                MeshMaterial3d(material_handle.clone()),
                Transform::default(),
                GlobalTransform::default(),
                BondVisualization {
                    atom1: bond.a,
                    atom2: bond.b,
                    strand_index: i,
                    total_strands: num_strands,
                },
            ));
        }
    }

    // Generate angles from the now-complete adjacency list
    for (center_atom, neighbors) in adjacency {
        if neighbors.len() < 2 {
            continue;
        }
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                connectivity.angles.push(Angle {
                    a: neighbors[i],
                    center: center_atom,
                    b: neighbors[j],
                });
            }
        }
    }

    info!(
        "Connectivity and visuals rebuilt. Bonds: {}, Angles: {}",
        connectivity.bonds.len(),
        connectivity.angles.len()
    );
}

#[derive(Resource, Default)]
pub struct SelectionState {
    pub selected: Vec<Entity>,
}

fn handle_selection(
    trigger: Trigger<Pointer<Click>>,
    mut selection: ResMut<SelectionState>,
    atoms: Query<&Atom>,
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
    // It also takes an EventWriter.
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
) {
    if keys.just_pressed(KeyCode::Delete) {
        let entities_to_delete = selection.selected.clone();
        if !entities_to_delete.is_empty() {
            info!(
                "Deletion key pressed for entities: {:?}",
                entities_to_delete
            );

            let mut connectivity_changed = false;
            for entity_to_delete in &entities_to_delete {
                commands.entity(*entity_to_delete).despawn();
                let initial_bond_count = connectivity.bonds.len();
                connectivity
                    .bonds
                    .retain(|bond| bond.a != *entity_to_delete && bond.b != *entity_to_delete);
                if connectivity.bonds.len() != initial_bond_count {
                    connectivity_changed = true;
                }
            }

            // Only trigger a rebuild if a bond was actually removed.
            if connectivity_changed {
                rebuild_writer.write(RebuildConnectivityEvent);
            }

            selection.selected.clear();
            info!("Cleared selection after deletion.");
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
        With<Atom>,
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
    atom_query: Query<(&Transform, &Force, &Velocity), With<Atom>>,
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
    atom_query: Query<&Transform, With<Atom>>,
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
