use bevy::pbr::wireframe::{Wireframe, WireframePlugin};
use bevy::prelude::*;
use bevy::ui::widget::TextUiWriter;
use bevy_picking::prelude::MeshPickingPlugin;
use bevy_picking::{
    events::{Click, Pointer},
    hover::PickingInteraction,
};

use crate::core::{AtomType, Force, Velocity};

use bevy::color::palettes::basic::{BLACK, BLUE, RED, YELLOW};

pub struct DebugPlugin;

impl Plugin for DebugPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(MeshPickingPlugin)
            .add_plugins(WireframePlugin::default())
            .init_resource::<SelectionState>()
            .add_systems(Startup, setup_debug_ui)
            .add_observer(handle_selection)
            .add_systems(
                Update,
                (update_highlights, update_info_panel, draw_selection_gizmos),
            );
    }
}

#[derive(Resource, Default)]
struct SelectionState {
    selected: Option<Entity>,
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
            let force_vector = force.0 * 0.0004;
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
    mut commands: Commands,
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
            let is_selected = Some(entity) == selection.selected;

            if is_selected {
                commands.entity(entity).insert(Wireframe);
                material.emissive = YELLOW.into();
            } else {
                material.emissive = BLACK.into();
                match interaction {
                    PickingInteraction::Hovered => {
                        commands.entity(entity).insert(Wireframe);
                    }
                    _ => {
                        commands.entity(entity).remove::<Wireframe>();
                    }
                }
            }
        }
    }
}

fn update_info_panel(
    selection: Res<SelectionState>,
    atom_query: Query<(&AtomType, &Transform, &Velocity, &Force)>,
    panel_query: Query<Entity, With<DebugInfoPanel>>,
    mut writer: TextUiWriter,
) {
    if let Ok(panel_entity) = panel_query.single() {
        if let Some(selected_entity) = selection.selected {
            if let Ok((atom_type, transform, velocity, force)) = atom_query.get(selected_entity) {
                *writer.text(panel_entity, 0) = format!(
                    "Selected Atom: {:?}\n\
                         Type: {:?}\n\
                         Position: {:.3?}\n\
                         Velocity (mag): {:.3?}\n\
                         Force (mag): {:.3?}",
                    selected_entity,
                    atom_type,
                    transform.translation,
                    velocity.length(),
                    force.length()
                );
            }
        } else {
            if selection.is_changed() {
                *writer.text(panel_entity, 0) = "Select an atom for details".to_string();
            }
        }
    }
}
