// src/ui/help_panel.rs

use crate::components::Atom;
use crate::interaction::SelectionState;
use bevy::picking::hover::PickingInteraction;
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

/// A resource to control the visibility of the help window.
#[derive(Resource, Default)]
struct HelpPanelState {
    is_open: bool,
}

pub struct HelpPanelPlugin;

impl Plugin for HelpPanelPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<HelpPanelState>()
            // We add a system to toggle visibility and one to draw the UI.
            .add_systems(Update, toggle_help_visibility)
            .add_systems(EguiPrimaryContextPass, help_panel_egui_system);
    }
}

/// This system listens for the 'H' key to toggle the help window's visibility.
fn toggle_help_visibility(keys: Res<ButtonInput<KeyCode>>, mut help_state: ResMut<HelpPanelState>) {
    if keys.just_pressed(KeyCode::KeyH) {
        help_state.is_open = !help_state.is_open;
    }
}

/// The single system that draws the entire help panel as an egui window.
fn help_panel_egui_system(
    mut contexts: EguiContexts,
    mut panel_state: ResMut<HelpPanelState>,
    selection: Res<SelectionState>,
    // We query for atoms being hovered over to provide contextual help.
    hover_query: Query<&PickingInteraction, With<Atom>>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    // This defines the window and binds its open/closed state to our resource.
    // The user can close it with the 'x' button, and 'H' will reopen it.
    if panel_state.is_open {
        egui::Window::new("Help")
            .open(&mut panel_state.is_open)
            .anchor(egui::Align2::RIGHT_BOTTOM, egui::vec2(-10.0, -10.0))
            .resizable(false)
            .show(ctx, |ui| {
                // Using a grid layout makes it easy to align key/description pairs.
                egui::Grid::new("help_grid")
                    .num_columns(2)
                    .spacing([20.0, 4.0])
                    .striped(true)
                    .show(ui, |ui| {
                        ui.label("h");
                        ui.label("Toggle this help window");
                        ui.end_row();

                        ui.label("Space");
                        ui.label("Toggle simulation pause/resume");
                        ui.end_row();

                        ui.label("LMB Drag");
                        ui.label("Rotate Camera");
                        ui.end_row();

                        ui.label("RMB Drag");
                        ui.label("Pan Camera");
                        ui.end_row();

                        ui.label("Scroll");
                        ui.label("Zoom Camera");
                        ui.end_row();

                        // Contextual help based on selection
                        if !selection.selected.is_empty() {
                            ui.separator();
                            ui.end_row();

                            ui.label("F");
                            ui.label("Focus camera on last selected atom");
                            ui.end_row();

                            ui.label("Shift+LMB");
                            ui.label("Add/Remove from multi-selection");
                            ui.end_row();

                            ui.label("Delete");
                            ui.label("Delete selected atom(s)");
                            ui.end_row();
                        }

                        if selection.selected.len() == 2 {
                            ui.separator();
                            ui.end_row();

                            ui.label("b");
                            ui.label("Create or delete bond");
                            ui.end_row();
                        }
                    });

                // Add a separator and some extra info at the bottom
                ui.separator();
                let is_hovering_atom = hover_query
                    .iter()
                    .any(|i| *i == PickingInteraction::Hovered);
                if is_hovering_atom {
                    ui.label("💡 Tip: Click on an atom to select it.");
                }
            });
    } else {
        let area_id = egui::Id::new("Help Opener Area");
        egui::Area::new(area_id)
            .anchor(egui::Align2::RIGHT_BOTTOM, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                // Add a frame for better visuals
                egui::Frame::popup(ui.style()).show(ui, |ui| {
                    // This button acts as the "stub" for the closed window.
                    let response = ui.button("❔ Help (H)");
                    if response.clicked() {
                        panel_state.is_open = true;
                    }
                    // Optional: Add a tooltip to the button itself
                    response.on_hover_text("Click or press 'H' to open the help panel.");
                });
            });
    }
}
