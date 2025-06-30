use crate::core::{AtomType, BondVisualization, Force, SystemConnectivity, Velocity};
use crate::simulation::PhysicsSet;
use bevy::color::palettes::basic::{BLACK, BLUE, RED, YELLOW};
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
    selected: Option<Entity>,
}

fn delete_selected_atom(
    mut commands: Commands,
    keys: Res<ButtonInput<KeyCode>>,
    mut selection: ResMut<SelectionState>,
    mut connectivity: ResMut<SystemConnectivity>,
    // We need to query for the bond visualization entities to despawn them too
    bond_vis_query: Query<(Entity, &BondVisualization)>,
) {
    // Check if the delete key was pressed and if an atom is selected
    if keys.just_pressed(KeyCode::Delete) {
        if let Some(entity_to_delete) = selection.selected {
            info!("Deletion key pressed for entity: {:?}", entity_to_delete);

            commands.entity(entity_to_delete).despawn();
            info!("Despawned atom entity.");

            connectivity.bonds.retain(|bond| {
                if bond.a == entity_to_delete || bond.b == entity_to_delete {
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

            // Keep only the angles that DO NOT contain the deleted entity in any position
            connectivity.angles.retain(|angle| {
                angle.a != entity_to_delete
                    && angle.b != entity_to_delete
                    && angle.center != entity_to_delete
            });

            info!("Cleaned up connectivity resource.");

            selection.selected = None;
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

    if let Some(selected_entity) = selection.selected {
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
) {
    if atoms.get(trigger.target()).is_ok() {
        if selection.selected == Some(trigger.target()) {
            selection.selected = None;
        } else {
            selection.selected = Some(trigger.target());
            info!("selected {:?}", selection.selected);
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

            if Some(entity) == selection.selected {
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
    atom_query: Query<(&AtomType, &Velocity, &Force)>,
    panel_query: Query<Entity, With<DebugInfoPanel>>,
    mut writer: TextUiWriter,
) {
    if let Ok(panel_entity) = panel_query.single() {
        if let Some(selected_entity) = selection.selected {
            if let Ok((atom_type, velocity, force)) = atom_query.get(selected_entity) {
                *writer.text(panel_entity, 0) = format!(
                    "Selected: {:?}\n\
                     Type: {:?}\n\
                     Vel (mag): {:.2}\n\
                     --- Forces (mag) ---\n\
                     Total: {:.1}\n\
                     Bond: {:.1}\n\
                     Angle: {:.1}\n\
                     Non-Bonded: {:.1}",
                    selected_entity,
                    atom_type,
                    velocity.length(),
                    force.total_magnitude(),
                    force.bond.length(),
                    force.angle.length(),
                    force.non_bonded.length()
                );
            }
        } else {
            if selection.is_changed() {
                *writer.text(panel_entity, 0) = "Select an atom for details".to_string();
            }
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
        if let Some(selected_entity) = selection.selected {
            if let Ok(mut pan_orbit_camera) = camera_query.single_mut() {
                if let Ok(atom_transform) = atom_query.get(selected_entity) {
                    pan_orbit_camera.target_focus = atom_transform.translation;
                }
            }
        }
    }
}
