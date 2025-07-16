// src/ui/pause_menu.rs

use crate::AppState;
use crate::components::Solvent;
use crate::persistence::{LoadStateEvent, SAVE_FILE_PATH, SaveStateEvent};
use crate::resources::{
    ForceMultiplier, LastSaveTime, PauseMenuState, SimulationParameters, Thermostat,
};
use crate::spawning::{DespawnSolventEvent, SpawnSolventEvent};
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
    mut next_state: ResMut<NextState<AppState>>,
    mut sim_params: ResMut<SimulationParameters>,
    mut thermostat: ResMut<Thermostat>,
    mut force_multiplier: ResMut<ForceMultiplier>,
    solvent_query: Query<Entity, With<Solvent>>,
    mut spawn_writer: EventWriter<SpawnSolventEvent>,
    mut despawn_writer: EventWriter<DespawnSolventEvent>,
    mut save_state_writer: EventWriter<SaveStateEvent>,
    mut load_state_writer: EventWriter<LoadStateEvent>,
    last_save_time: Res<LastSaveTime>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::Area::new(egui::Id::new("settings_button_area"))
        .anchor(egui::Align2::RIGHT_TOP, egui::vec2(-10.0, 10.0))
        .show(ctx, |ui| {
            if !menu_state.visible {
                ui.horizontal(|ui| {
                    if ui.button("Settings (M)").clicked() {
                        menu_state.visible = true;
                        next_state.set(AppState::Paused);
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

            // --- ADD THIS LOGIC ---
            let is_solvent_present = !solvent_query.is_empty();
            let button_text = if is_solvent_present {
                "Despawn Solvent"
            } else {
                "Spawn Solvent"
            };

            if ui.button(button_text).clicked() {
                if is_solvent_present {
                    despawn_writer.write(DespawnSolventEvent);
                } else {
                    spawn_writer.write(SpawnSolventEvent);
                }
            }
            #[cfg(not(target_arch = "wasm32"))]
            {
                ui.add_space(10.0);
                ui.separator();
                ui.add_space(10.0);

                ui.horizontal(|ui| {
                    if ui.button("Save State").clicked() {
                        save_state_writer.write(SaveStateEvent(SAVE_FILE_PATH.to_string()));
                    }
                    if ui.button("Load State").clicked() {
                        load_state_writer.write(LoadStateEvent(SAVE_FILE_PATH.to_string()));
                    }
                });
                ui.add_space(2.0);
                ui.label(&last_save_time.display_text);
            }
            ui.add_space(10.0);
            ui.separator();

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
