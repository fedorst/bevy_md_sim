// src/cursor.rs

use crate::components::Atom;
use crate::interaction::SelectionState;
use crate::resources::DragState;
use bevy::prelude::*;
use bevy::winit::cursor::CursorIcon;
use bevy_egui::EguiContexts;
use bevy_picking::hover::PickingInteraction;

// An enum to represent the logical state of the cursor.
#[derive(Debug, PartialEq, Eq, Clone, Copy, Default, States, Hash)]
enum CursorState {
    #[default]
    Default,
    Hover,
    Dragging,
    Draggable,
}

pub struct CustomCursorPlugin;

impl Plugin for CustomCursorPlugin {
    fn build(&self, app: &mut App) {
        app.init_state::<CursorState>()
            .add_systems(Startup, setup_cursor_icons)
            .add_systems(
                Update,
                (
                    determine_cursor_state,
                    apply_cursor_icon.run_if(state_changed::<CursorState>),
                )
                    .chain(),
            );
    }
}

// Resource to store our cursor icon handles
#[derive(Resource)]
struct CursorIcons {
    default: CursorIcon,
    hover: CursorIcon,
    dragging: CursorIcon,
    draggable: CursorIcon,
}

fn setup_cursor_icons(mut commands: Commands, asset_server: Res<AssetServer>) {
    // All custom cursors using images
    let default_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle: asset_server.load("cursors/pointer_b.png"),
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: (16, 16), // Adjust based on your image
        },
    ));

    let hover_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle: asset_server.load("cursors/pointer_c.png"),
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: (16, 16), // Adjust based on your image
        },
    ));

    let draggable_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle: asset_server.load("cursors/hand_open.png"),
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: (16, 16), // Adjust based on your image
        },
    ));

    let dragging_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle: asset_server.load("cursors/hand_closed.png"),
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: (16, 16), // Adjust based on your image
        },
    ));

    commands.insert_resource(CursorIcons {
        default: default_cursor,
        hover: hover_cursor,
        dragging: dragging_cursor,
        draggable: draggable_cursor,
    });
}

/// This system determines what the cursor's state SHOULD be.
fn determine_cursor_state(
    drag_state: Res<DragState>,
    mut egui_contexts: EguiContexts,
    atom_hover_q: Query<(Entity, &PickingInteraction), With<Atom>>,
    selection: Res<SelectionState>,
    mut next_state: ResMut<NextState<CursorState>>,
) {
    let Ok(ctx) = egui_contexts.ctx_mut() else {
        return;
    };

    let mut new_state = CursorState::Default;

    // Priority 1: A drag is officially in progress (mouse has moved).
    if !drag_state.dragged_entities.is_empty() {
        new_state = CursorState::Dragging;
    }
    // Priority 2: The cursor is over an egui element.
    else if ctx.is_pointer_over_area() || ctx.is_using_pointer() {
        new_state = CursorState::Hover;
    }
    // Priority 3: Check for interactions in the 3D scene.
    else {
        let mut is_hovering_any_atom = false;
        let mut is_hovering_selected_atom = false;
        let mut is_pressing_selected_atom = false;

        for (entity, interaction) in &atom_hover_q {
            if *interaction != PickingInteraction::None {
                is_hovering_any_atom = true;
                // Is the interacted atom part of our selection?
                if !selection.selected.is_empty() && selection.selected.contains(&entity) {
                    match *interaction {
                        PickingInteraction::Pressed => {
                            is_pressing_selected_atom = true;
                            // This is the highest priority scene interaction, so we can stop looking.
                            break;
                        }
                        PickingInteraction::Hovered => {
                            is_hovering_selected_atom = true;
                        }
                        _ => {}
                    }
                }
            }
        }

        // Now, set the state based on the flags we collected, in order of priority.
        if is_pressing_selected_atom {
            new_state = CursorState::Dragging; // A press on a selected item IS dragging.
        } else if is_hovering_selected_atom {
            new_state = CursorState::Draggable; // A hover on a selected item means it's ready to drag.
        } else if is_hovering_any_atom {
            new_state = CursorState::Hover; // A hover on any other item is just a hover.
        }
        // If no flags are set, the state remains CursorState::Default.
    }

    next_state.set(new_state);
}

/// This system runs ONLY when the CursorState changes and applies the new cursor icon.
fn apply_cursor_icon(
    mut windows: Query<Entity, With<Window>>,
    state: Res<State<CursorState>>,
    cursor_icons: Res<CursorIcons>,
    mut commands: Commands,
) {
    let Ok(window_entity) = windows.single_mut() else {
        return;
    };

    let new_cursor = match state.get() {
        CursorState::Default => &cursor_icons.default,
        CursorState::Hover => &cursor_icons.hover,
        CursorState::Dragging => &cursor_icons.dragging,
        CursorState::Draggable => &cursor_icons.draggable,
    };

    commands.entity(window_entity).insert(new_cursor.clone());
}
