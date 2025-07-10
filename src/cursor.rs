// src/cursor.rs

use crate::components::Atom;
use crate::cursor;
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
            .add_systems(Startup, load_cursor_assets)
            .add_systems(
                Update,
                (
                    create_resized_cursors
                        .run_if(resource_exists_and_changed::<CursorAssetHandles>),
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

#[derive(Resource)]
struct CursorAssetHandles {
    default: Handle<Image>,
    hover: Handle<Image>,
    dragging: Handle<Image>,
    draggable: Handle<Image>,
}

fn load_cursor_assets(mut commands: Commands, asset_server: Res<AssetServer>) {
    commands.insert_resource(CursorAssetHandles {
        default: asset_server.load("cursors/pointer_b.png"),
        hover: asset_server.load("cursors/pointer_c.png"),
        dragging: asset_server.load("cursors/hand_closed.png"),
        draggable: asset_server.load("cursors/hand_open.png"),
    });
}

fn create_resized_cursors(
    mut commands: Commands,
    handles: Res<CursorAssetHandles>,
    mut images: ResMut<Assets<Image>>,
) {
    if let (Some(default_img), Some(hover_img), Some(dragging_img), Some(draggable_img)) = (
        images.get(&handles.default),
        images.get(&handles.hover),
        images.get(&handles.dragging),
        images.get(&handles.draggable),
    ) {
        // Resize images and create new handles
        let default_resized = resize_image(default_img, 32, 32);
        let hover_resized = resize_image(hover_img, 32, 32);
        let dragging_resized = resize_image(dragging_img, 32, 32);
        let draggable_resized = resize_image(draggable_img, 32, 32);

        let default_handle = images.add(default_resized);
        let hover_handle = images.add(hover_resized);
        let dragging_handle = images.add(dragging_resized);
        let draggable_handle = images.add(draggable_resized);

        // Create cursor icons
        let default_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
            bevy::winit::cursor::CustomCursorImage {
                handle: default_handle,
                texture_atlas: None,
                flip_x: false,
                flip_y: false,
                rect: None,
                hotspot: (0, 0),
            },
        ));

        let hover_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
            bevy::winit::cursor::CustomCursorImage {
                handle: hover_handle,
                texture_atlas: None,
                flip_x: false,
                flip_y: false,
                rect: None,
                hotspot: (0, 0),
            },
        ));

        let dragging_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
            bevy::winit::cursor::CustomCursorImage {
                handle: dragging_handle,
                texture_atlas: None,
                flip_x: false,
                flip_y: false,
                rect: None,
                hotspot: (0, 0),
            },
        ));

        let draggable_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
            bevy::winit::cursor::CustomCursorImage {
                handle: draggable_handle,
                texture_atlas: None,
                flip_x: false,
                flip_y: false,
                rect: None,
                hotspot: (0, 0),
            },
        ));

        commands.insert_resource(CursorIcons {
            default: default_cursor,
            hover: hover_cursor,
            dragging: dragging_cursor,
            draggable: draggable_cursor,
        });

        commands.remove_resource::<CursorAssetHandles>();
        info!("Resized cursor images created");
    }
}

fn resize_image(image: &Image, new_width: u32, new_height: u32) -> Image {
    let old_width = image.texture_descriptor.size.width;
    let old_height = image.texture_descriptor.size.height;

    if old_width == new_width && old_height == new_height {
        return image.clone();
    }

    let mut new_data = Vec::with_capacity((new_width * new_height * 4) as usize);

    let Some(old_data) = &image.data else {
        return image.clone();
    };

    for y in 0..new_height {
        for x in 0..new_width {
            // Nearest neighbor sampling
            let src_x = (x * old_width / new_width).min(old_width - 1);
            let src_y = (y * old_height / new_height).min(old_height - 1);
            let src_index = ((src_y * old_width + src_x) * 4) as usize;

            // Copy RGBA pixels
            if src_index + 3 < old_data.len() {
                new_data.extend_from_slice(&old_data[src_index..src_index + 4]);
            } else {
                new_data.extend_from_slice(&[0, 0, 0, 0]); // Transparent fallback
            }
        }
    }

    let mut new_image = image.clone();
    new_image.data = Some(new_data);
    new_image.texture_descriptor.size.width = new_width;
    new_image.texture_descriptor.size.height = new_height;

    new_image
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
            hotspot: (0, 0), // Adjust based on your image
        },
    ));

    let hover_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle: asset_server.load("cursors/pointer_c.png"),
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: (0, 0), // Adjust based on your image
        },
    ));

    let draggable_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle: asset_server.load("cursors/hand_open.png"),
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: (0, 0), // Adjust based on your image
        },
    ));

    let dragging_cursor = CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle: asset_server.load("cursors/hand_closed.png"),
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: (0, 0), // Adjust based on your image
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
