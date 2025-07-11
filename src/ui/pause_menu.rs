// src/ui/pause_menu.rs

use crate::resources::{
    ForceMultiplier, PauseMenuState, SimulationParameters, SimulationState, Thermostat,
};
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

pub struct PauseMenuPlugin;

impl Plugin for PauseMenuPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Update, toggle_pause_menu_visibility)
            .add_systems(EguiPrimaryContextPass, pause_menu_egui_system);
    }
}

fn toggle_pause_menu_visibility(
    keys: Res<ButtonInput<KeyCode>>,
    mut menu_state: ResMut<PauseMenuState>,
) {
    if keys.just_pressed(KeyCode::KeyM) {
        menu_state.visible = !menu_state.visible;
    }
}

/// The single system that draws the entire pause menu as a modal egui window.
fn pause_menu_egui_system(
    mut contexts: EguiContexts,
    mut menu_state: ResMut<PauseMenuState>,
    mut sim_state: ResMut<SimulationState>,
    mut sim_params: ResMut<SimulationParameters>,
    mut thermostat: ResMut<Thermostat>,
    mut force_multiplier: ResMut<ForceMultiplier>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::Area::new(egui::Id::new("settings_button_area"))
        .anchor(egui::Align2::RIGHT_TOP, egui::vec2(-10.0, 10.0))
        .show(ctx, |ui| {
            // If the menu is closed, show the button to open it.
            if !menu_state.visible {
                // THE FIX for the layout bug: Use a horizontal layout to prevent wrapping.
                ui.horizontal(|ui| {
                    if ui.button("Settings (M)").clicked() {
                        menu_state.visible = true;
                        sim_state.paused = true;
                    }
                });
            }
        });

    if !menu_state.visible {
        return;
    }
    let mut is_open = menu_state.visible;
    egui::Window::new("Settings")
        .open(&mut is_open) // Binds to the 'x' button
        .collapsible(false)
        .resizable(false) // Keep it a fixed size
        .default_pos(ctx.screen_rect().center() - egui::vec2(150.0, 300.0)) // Approx center
        .id(egui::Id::new("settings_window"))
        .show(ctx, |ui| {
            ui.set_width(300.0);

            egui::Grid::new("settings_grid")
                .num_columns(2)
                .spacing([40.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Timestep (dt)");
                    ui.add(
                        egui::DragValue::new(&mut sim_params.dt)
                            .speed(0.00001)
                            .range(0.0..=0.1)
                            .prefix("ps: "),
                    );
                    ui.end_row();

                    ui.label("Target Temp");
                    ui.add(
                        egui::DragValue::new(&mut thermostat.target_temperature)
                            .speed(1.0)
                            .range(0.0..=1000.0)
                            .suffix(" K"),
                    );
                    ui.end_row();

                    ui.label("Thermostat Tau");
                    ui.add(
                        egui::DragValue::new(&mut thermostat.tau)
                            .speed(0.001)
                            .range(0.0..=1.0)
                            .prefix("ps: "),
                    );
                    ui.end_row();

                    ui.label("Force Multiplier");
                    ui.add(
                        egui::DragValue::new(&mut force_multiplier.0)
                            .speed(0.05)
                            .range(0.0..=10.0),
                    );
                    ui.end_row();
                });

            ui.add_space(10.0);
            ui.separator();
            ui.add_space(10.0);

            ui.with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
                if ui.button("Close").clicked() {
                    menu_state.visible = false;
                }
            });
        });

    // If the user closed the window with the 'x' button, update our state.
    if !is_open {
        menu_state.visible = false;
    }
}
