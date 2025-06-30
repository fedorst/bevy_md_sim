use crate::core::{
    AtomIdMap, AtomType, BondVisualization, Force, ForceField, SystemConnectivity, Velocity,
};
use crate::simulation::PhysicsSet;
use bevy::color::palettes::basic::{BLACK, BLUE, RED};
use bevy::prelude::*;
use bevy::ui::widget::TextUiWriter;
use bevy_panorbit_camera::PanOrbitCamera;
use bevy_picking::prelude::MeshPickingPlugin;
use bevy_picking::{
    events::{Click, Pointer},
    hover::PickingInteraction,
};
pub struct DebugPlugin;

impl Plugin for DebugPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(MeshPickingPlugin)
            .init_resource::<SelectionState>()
            .add_systems(Startup, setup_debug_ui)
            .add_observer(handle_selection)
            .add_systems(
                Update,
                (
                    update_highlights,
                    update_info_panel,
                    draw_selection_gizmos,
                    delete_selected_atom,
                    focus_camera_on_selection,
                )
                    .after(PhysicsSet),
            );
    }
}

#[derive(Resource, Default)]
struct SelectionState {
    selected: Vec<Entity>,
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

#[derive(Component)]
struct DebugInfoPanel;

fn setup_debug_ui(mut commands: Commands) {
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            bottom: Val::Px(10.0),
            left: Val::Px(10.0),
            width: Val::Px(250.0),
            padding: UiRect::all(Val::Px(10.0)),
            flex_direction: FlexDirection::Column,
            ..default()
        },
        BackgroundColor(Color::BLACK.with_alpha(0.75)),
        Text::new("Select an atom for details"),
        TextFont {
            font_size: 16.0,
            ..default()
        },
        DebugInfoPanel,
    ));
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

fn update_info_panel(
    selection: Res<SelectionState>,
    connectivity: Res<SystemConnectivity>,
    force_field: Res<ForceField>,
    atom_id_map: Res<AtomIdMap>,
    atom_query: Query<(&AtomType, &Transform, &Velocity, &Force)>,
    panel_query: Query<Entity, With<DebugInfoPanel>>,
    mut writer: TextUiWriter,
) {
    if let Ok(panel_entity) = panel_query.single() {
        let mut info_text = "Select an atom for details".to_string();

        match selection.selected.len() {
            0 => {
                // Text is already set to default.
            }
            1 => {
                // --- SINGLE ATOM DISPLAY (as before) ---
                let entity = selection.selected[0];
                if let Ok((atom_type, _transform, velocity, force)) = atom_query.get(entity) {
                    let pretty_id = atom_id_map.entity_to_id.get(&entity).unwrap();
                    info_text = format!(
                        "Selected: {}\n\
                         Type: {:?}\n\
                         Vel (mag): {:.2}\n\
                         --- Forces (mag) ---\n\
                         Total: {:.1}\n\
                         Bond: {:.1}\n\
                         Angle: {:.1}\n\
                         Non-Bonded: {:.1}",
                        pretty_id,
                        atom_type,
                        velocity.length(),
                        force.total_magnitude(),
                        force.bond.length(),
                        force.angle.length(),
                        force.non_bonded.length()
                    );
                }
            }
            2 => {
                let e1 = selection.selected[0];
                let e2 = selection.selected[1];

                if let Some(_) = connectivity
                    .bonds
                    .iter()
                    .find(|b| (b.a == e1 && b.b == e2) || (b.a == e2 && b.b == e1))
                {
                    if let Ok([(type1, t1, _, _), (type2, t2, _, _)]) =
                        atom_query.get_many([e1, e2])
                    {
                        if let Some(&(_k, r0)) = force_field.bond_params.get(&(*type1, *type2)) {
                            let current_dist = t1.translation.distance(t2.translation);
                            let id1 = atom_id_map.entity_to_id.get(&e1).unwrap();
                            let id2 = atom_id_map.entity_to_id.get(&e2).unwrap();
                            info_text = format!(
                                "Selected Bond:\n{}-{}\n\
                                 --- Geometry ---\n\
                                 Current Len: {:.4} nm\n\
                                 Optimal Len: {:.4} nm\n\
                                 Strain: {:.2}%",
                                id1,
                                id2,
                                current_dist,
                                r0,
                                (current_dist / r0 - 1.0) * 100.0 // Strain as a percentage
                            );
                        }
                    }
                } else {
                    // Default multi-select text if no bond.
                    info_text = format!("Selected {} atoms.", selection.selected.len());
                }
            }
            3 => {
                let entities = &selection.selected;
                if let Some(angle) = connectivity.angles.iter().find(|a| {
                    entities.contains(&a.a)
                        && entities.contains(&a.center)
                        && entities.contains(&a.b)
                }) {
                    // Get the data for all three atoms
                    if let Ok(
                        [
                            (type_a, t_a, _, _),
                            (type_center, t_center, _, _),
                            (type_b, t_b, _, _),
                        ],
                    ) = atom_query.get_many([angle.a, angle.center, angle.b])
                    {
                        // Get the ideal angle (theta0) from the force field
                        if let Some(&(_k, theta0)) =
                            force_field
                                .angle_params
                                .get(&(*type_a, *type_center, *type_b))
                        {
                            let v1 = t_a.translation - t_center.translation;
                            let v2 = t_b.translation - t_center.translation;
                            let current_angle_rad = v1.angle_between(v2);
                            let id_a = atom_id_map.entity_to_id.get(&angle.a).unwrap();
                            let id_center = atom_id_map.entity_to_id.get(&angle.center).unwrap();
                            let id_b = atom_id_map.entity_to_id.get(&angle.b).unwrap();
                            info_text = format!(
                                "Selected Angle:\n{}-{}-{}\n\
                                 --- Geometry ---\n\
                                 Optimal: {:.2}deg\n\
                                 Current: {:.2}deg ({:+.2})\n",
                                id_a,
                                id_center,
                                id_b,
                                theta0.to_degrees(),
                                current_angle_rad.to_degrees(),
                                current_angle_rad.to_degrees() - theta0.to_degrees()
                            );
                        }
                    }
                } else {
                    info_text = format!("Selected {} atoms.", selection.selected.len());
                }
            }
            _ => {
                info_text = format!("Selected {} atoms.", selection.selected.len());
            }
        }
        *writer.text(panel_entity, 0) = info_text;
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
