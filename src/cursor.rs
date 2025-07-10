// src/cursor.rs

use crate::components::Atom;
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
    atom_hover_q: Query<&PickingInteraction, With<Atom>>,
    mut next_state: ResMut<NextState<CursorState>>,
) {
    let Ok(ctx) = egui_contexts.ctx_mut() else {
        return;
    };

    let mut new_state = CursorState::Default;

    if !drag_state.dragged_entities.is_empty() {
        new_state = CursorState::Dragging;
    } else if ctx.is_pointer_over_area() || ctx.is_using_pointer() {
        new_state = CursorState::Hover;
    } else if atom_hover_q.iter().any(|i| *i != PickingInteraction::None) {
        new_state = CursorState::Hover;
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
    };

    commands.entity(window_entity).insert(new_cursor.clone());
}
