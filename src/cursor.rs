// src/cursor.rs

use crate::components::Atom;
use crate::interaction::SelectionState;
use crate::resources::DragState;
use bevy::input::mouse::MouseWheel;
use bevy::prelude::*;
use bevy::winit::cursor::CursorIcon;
use bevy_egui::EguiContexts;
use bevy_picking::hover::PickingInteraction;

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq, Hash, States)]
enum CursorLoadingState {
    #[default]
    NotLoaded,
    Loading,
    Finished,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Default, States, Hash)]
enum CursorState {
    #[default]
    Default,
    Hover,
    Dragging,
    Draggable,
    ZoomingIn,
    ZoomingOut,
    // the rest is not in use until i figure out how to override egui pointer preferences
    Text,
    ResizeEW, // East-West
    ResizeNS, // North-South
    ResizeNESW,
    ResizeNWSE,
}

pub struct CustomCursorPlugin;

impl Plugin for CustomCursorPlugin {
    fn build(&self, app: &mut App) {
        app.init_state::<CursorLoadingState>()
            .init_state::<CursorState>()
            .init_resource::<ZoomStateTimer>()
            .add_systems(Startup, load_cursor_assets)
            .add_systems(
                Update,
                create_resized_cursors.run_if(in_state(CursorLoadingState::Loading)),
            )
            .add_systems(
                Update,
                (
                    handle_zoom_cursor,
                    determine_cursor_state,
                    apply_cursor_icon.run_if(in_state(CursorLoadingState::Finished)),
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
    text: CursorIcon,
    zoom_in: CursorIcon,
    zoom_out: CursorIcon,
    resize_ew: CursorIcon,
    resize_ns: CursorIcon,
    resize_nesw: CursorIcon,
    resize_nwse: CursorIcon,
}

#[derive(Resource)]
struct CursorAssetHandles {
    default: Handle<Image>,
    hover: Handle<Image>,
    dragging: Handle<Image>,
    draggable: Handle<Image>,
    text: Handle<Image>,
    zoom_in: Handle<Image>,
    zoom_out: Handle<Image>,
    resize_ew: Handle<Image>,
    resize_ns: Handle<Image>,
    resize_nesw: Handle<Image>,
    resize_nwse: Handle<Image>,
}

#[derive(Resource)]
struct ZoomStateTimer(Timer);

impl Default for ZoomStateTimer {
    fn default() -> Self {
        let mut timer = Timer::from_seconds(0.2, TimerMode::Once);
        timer.pause();
        Self(timer)
    }
}

fn load_cursor_assets(
    mut commands: Commands,
    asset_server: Res<AssetServer>,
    mut next_state: ResMut<NextState<CursorLoadingState>>,
) {
    info!("[Cursor] Starting to load cursor image assets...");
    commands.insert_resource(CursorAssetHandles {
        default: asset_server.load("cursors/pointer_b.png"),
        hover: asset_server.load("cursors/pointer_c.png"),
        dragging: asset_server.load("cursors/hand_closed.png"),
        draggable: asset_server.load("cursors/hand_open.png"),
        text: asset_server.load("cursors/bracket_a_vertical.png"),
        resize_ew: asset_server.load("cursors/resize_a_horizontal.png"),
        resize_ns: asset_server.load("cursors/resize_a_vertical.png"),
        resize_nesw: asset_server.load("cursors/resize_a_diagonal.png"),
        resize_nwse: asset_server.load("cursors/resize_a_diagonal_mirror.png"),
        zoom_in: asset_server.load("cursors/zoom_in.png"),
        zoom_out: asset_server.load("cursors/zoom_out.png"),
    });
    next_state.set(CursorLoadingState::Loading);
}

fn process_cursor_image(
    image: &Image,
    hotspot_xy: (u16, u16),
    images: &mut ResMut<Assets<Image>>,
) -> CursorIcon {
    // 1. Resize the image
    let resized_image = resize_image(image, 32, 32);
    // 2. Add the resized image to the asset collection and get its new handle
    let handle = images.add(resized_image);

    // 3. Create and return the final CursorIcon
    CursorIcon::Custom(bevy::winit::cursor::CustomCursor::Image(
        bevy::winit::cursor::CustomCursorImage {
            handle,
            texture_atlas: None,
            flip_x: false,
            flip_y: false,
            rect: None,
            hotspot: hotspot_xy,
        },
    ))
}

fn create_resized_cursors(
    mut commands: Commands,
    handles: Res<CursorAssetHandles>,
    mut images: ResMut<Assets<Image>>,
    mut next_state: ResMut<NextState<CursorLoadingState>>,
) {
    // waiting until https://github.com/bevyengine/bevy/issues/17276 is done, for now manual resize

    // Define all our cursor data in one clean place.
    // Format: (The handle from CursorAssetHandles, (hotspot_x, hotspot_y))
    let cursor_definitions = [
        (&handles.default, (0, 0)),
        (&handles.hover, (0, 0)),
        (&handles.dragging, (0, 0)),  // Center the hotspot for grabbing
        (&handles.draggable, (0, 0)), // Center the hotspot for grabbing
        (&handles.text, (0, 0)),
        (&handles.zoom_in, (8, 8)),
        (&handles.zoom_out, (8, 8)),
        (&handles.resize_ew, (0, 0)),
        (&handles.resize_ns, (0, 0)),
        (&handles.resize_nesw, (0, 0)),
        (&handles.resize_nwse, (0, 0)),
    ];

    // First, check if all assets are loaded (this logic is from the state-based fix)
    if cursor_definitions
        .iter()
        .all(|(handle, _)| images.get(*handle).is_some())
    {
        info!("[Cursor] All cursor assets are loaded. Processing now.");

        // Since we know all images exist, we can now safely process them in a loop.
        let loaded_images: Vec<Image> = cursor_definitions
            .iter()
            .map(|(handle, _)| images.get(*handle).unwrap().clone())
            .collect();
        let processed_cursors: Vec<CursorIcon> = loaded_images
            .iter()
            .zip(cursor_definitions.iter())
            .map(|(image, (_, hotspot))| process_cursor_image(image, *hotspot, &mut images))
            .collect();

        commands.insert_resource(CursorIcons {
            default: processed_cursors[0].clone(),
            hover: processed_cursors[1].clone(),
            dragging: processed_cursors[2].clone(),
            draggable: processed_cursors[3].clone(),
            text: processed_cursors[4].clone(),
            zoom_in: processed_cursors[5].clone(),
            zoom_out: processed_cursors[6].clone(),
            resize_ew: processed_cursors[7].clone(),
            resize_ns: processed_cursors[8].clone(),
            resize_nesw: processed_cursors[9].clone(),
            resize_nwse: processed_cursors[10].clone(),
        });

        // Transition to the Finished state and clean up.
        next_state.set(CursorLoadingState::Finished);
        commands.remove_resource::<CursorAssetHandles>();
        info!("Resized cursor images created and resource inserted.");
    }
}

fn handle_zoom_cursor(
    mut scroll_evr: EventReader<MouseWheel>,
    mut next_state: ResMut<NextState<CursorState>>,
    mut zoom_timer: ResMut<ZoomStateTimer>,
) {
    if !scroll_evr.is_empty() {
        // Since we know we are scrolling, unpause and reset the timer now.
        zoom_timer.0.unpause();
        zoom_timer.0.reset();

        // Use a fold to get the net scroll direction for the frame
        let net_y = scroll_evr.read().fold(0.0, |acc, ev| acc + ev.y);

        if net_y > 0.0 {
            next_state.set(CursorState::ZoomingIn);
        } else if net_y < 0.0 {
            next_state.set(CursorState::ZoomingOut);
        }
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

/// This system determines what the cursor's state SHOULD be.
fn determine_cursor_state(
    drag_state: Res<DragState>,
    mut egui_contexts: EguiContexts,
    atom_hover_q: Query<(Entity, &PickingInteraction), With<Atom>>,
    selection: Res<SelectionState>,
    mut next_state: ResMut<NextState<CursorState>>,
    mut zoom_timer: ResMut<ZoomStateTimer>,
    time: Res<Time>,
    current_cursor_state: Res<State<CursorState>>,
) {
    let is_zooming = matches!(
        current_cursor_state.get(),
        CursorState::ZoomingIn | CursorState::ZoomingOut
    );

    if is_zooming {
        zoom_timer.0.tick(time.delta());
        // If the timer finishes, it's time to go back to default.
        if zoom_timer.0.finished() {
            info!("Finished zooming!");
            next_state.set(CursorState::Default);
        }
        // If we are zooming, we don't care about any other cursor logic.
        return;
    }

    let Ok(ctx) = egui_contexts.ctx_mut() else {
        info!("Egui not yet initialized for cursor determination");
        return;
    };

    let mut new_state = CursorState::Default;

    // Priority 1: A drag is officially in progress (mouse has moved).
    if !drag_state.dragged_entities.is_empty() {
        new_state = CursorState::Dragging;
    }
    // Priority 2: The cursor is over an egui element.
    else if ctx.is_pointer_over_area() || ctx.is_using_pointer() {
        return;
        // new_state = CursorState::Hover;
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
    current_cursor_state: Res<State<CursorState>>,
    cursor_icons: Res<CursorIcons>,
    mut commands: Commands,
    mut last_applied_state: Local<Option<CursorState>>,
) {
    let Ok(window_entity) = windows.single_mut() else {
        return;
    };
    let target_state = *current_cursor_state.get();

    if last_applied_state.is_none() || *last_applied_state.as_ref().unwrap() != target_state {
        let new_cursor_icon = match target_state {
            CursorState::Default => &cursor_icons.default,
            CursorState::Hover => &cursor_icons.hover,
            CursorState::Dragging => &cursor_icons.dragging,
            CursorState::Draggable => &cursor_icons.draggable,
            CursorState::Text => &cursor_icons.text,
            CursorState::ZoomingIn => &cursor_icons.zoom_in,
            CursorState::ZoomingOut => &cursor_icons.zoom_out,
            CursorState::ResizeEW => &cursor_icons.resize_ew,
            CursorState::ResizeNS => &cursor_icons.resize_ns,
            CursorState::ResizeNESW => &cursor_icons.resize_nesw,
            CursorState::ResizeNWSE => &cursor_icons.resize_nwse,
        };

        commands
            .entity(window_entity)
            .insert(new_cursor_icon.clone());
        *last_applied_state = Some(target_state);
        info!("CURSOR: Applied state {:?}", target_state); // Optional debug log
    }
}
