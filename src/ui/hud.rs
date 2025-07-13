// src/ui/hud.rs

use crate::AppState;
use crate::resources::{
    ActiveWallTime, CurrentTemperature, SimulationParameters, StepCount, SystemEnergy, Thermostat,
    ThermostatScale,
};
use crate::spawning::{SMILESValidationResult, SpawnMoleculeFromSMILESEvent, ValidateSMILESEvent};
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

#[derive(Resource, Default)]
struct EguiInputState {
    smiles: String,
    is_valid: bool,
    error_message: String,
}

pub struct HudPlugin;

impl Plugin for HudPlugin {
    fn build(&self, app: &mut App) {
        // The plugin now only needs to register one system.
        app.init_resource::<EguiInputState>()
            .add_systems(EguiPrimaryContextPass, hud_egui_system)
            .add_systems(Update, update_validation_state);
    }
}

/// This single system draws all the heads-up display elements using egui.
fn hud_egui_system(
    mut contexts: EguiContexts,
    // Query for all the data resources needed for the display.
    mut egui_state: ResMut<EguiInputState>,
    mut validate_writer: EventWriter<ValidateSMILESEvent>,
    mut spawn_writer: EventWriter<SpawnMoleculeFromSMILESEvent>,
    app_state: Res<State<AppState>>,
    active_wall_time: Res<ActiveWallTime>,
    step_count: Res<StepCount>,
    sim_params: Res<SimulationParameters>,
    energy: Res<SystemEnergy>,
    current_temp: Res<CurrentTemperature>,
    thermostat: Res<Thermostat>,
    thermostat_scale: Res<ThermostatScale>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    // --- 1. Top-Left Info Panel (Time, Energy, Temp) ---
    egui::Area::new(egui::Id::new("hud_info_area"))
        .anchor(egui::Align2::LEFT_TOP, egui::vec2(10.0, 10.0))
        .show(ctx, |ui| {
            // Use a semitransparent frame for readability
            let frame = egui::Frame::popup(ui.style()).fill(egui::Color32::from_black_alpha(128));
            frame.show(ui, |ui| {
                ui.set_width(220.0); // Give it a fixed width
                ui.label(
                    egui::RichText::new("Simulation Info").font(egui::FontId::proportional(16.0)),
                );
                ui.separator();
                let simulated_time_ps = step_count.0 as f32 * sim_params.dt;
                ui.label(format!("Sim Time: {:.3} ps", simulated_time_ps));
                ui.label(format!("Wall Time: {:.2} s", active_wall_time.0));
                ui.label(format!("Steps: {}", step_count.0));
                ui.separator();
                ui.label(format!(
                    "Temp: {:.1} K / {:.1} K",
                    current_temp.0, thermostat.target_temperature
                ));
                ui.label(format!("Scale: {:.6}", thermostat_scale.0));
                ui.separator();
                ui.label(format!("Total E: {:.2}", energy.total));
                ui.label(format!("Potential E: {:.2}", energy.potential));
                ui.label(format!("Kinetic E: {:.2}", energy.kinetic));
            });
        });

    // --- 2. Center "PAUSED" / "RUNNING" Text ---
    // Only draw this if the simulation is paused.
    if *app_state.get() == AppState::Paused {
        egui::Area::new(egui::Id::new("hud_paused_text"))
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, -100.0))
            .show(ctx, |ui| {
                let text = egui::RichText::new("PAUSED")
                    .font(egui::FontId::proportional(48.0))
                    .color(egui::Color32::from_white_alpha(180))
                    .strong();
                ui.label(text);
            });
    }

    let smiles_window = egui::Window::new("Molecule Input")
        .anchor(egui::Align2::CENTER_BOTTOM, egui::vec2(0.0, -20.0))
        .collapsible(false)
        .resizable(false)
        .title_bar(false);

    let frame_stroke = if egui_state.is_valid {
        egui::Stroke::new(1.0, egui::Color32::from_rgb(0, 255, 255)) // Cyan
    } else {
        egui::Stroke::new(2.0, egui::Color32::from_rgb(255, 50, 50)) // Red
    };
    let frame = egui::Frame::group(&ctx.style()).stroke(frame_stroke);

    smiles_window.frame(frame).show(ctx, |ui| {
        ui.set_width(300.0);

        let text_edit = egui::TextEdit::singleline(&mut egui_state.smiles)
            .hint_text("Type SMILES and press Enter...");

        let mut response = ui.add(text_edit);

        if !egui_state.is_valid {
            response = response.on_hover_text(&egui_state.error_message);
        }

        if response.changed() {
            validate_writer.write(ValidateSMILESEvent(egui_state.smiles.clone()));
        }

        if response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
            if !egui_state.smiles.is_empty() {
                spawn_writer.write(SpawnMoleculeFromSMILESEvent(
                    egui_state.smiles.clone(),
                    None,
                ));
            }
        }
    });
}

fn update_validation_state(
    mut egui_state: ResMut<EguiInputState>,
    mut validation_results: EventReader<SMILESValidationResult>,
) {
    if let Some(event) = validation_results.read().last() {
        match &event.0 {
            Ok(_) => {
                egui_state.is_valid = true;
                egui_state.error_message.clear();
            }
            Err(e) => {
                egui_state.is_valid = false;
                if e.trim().is_empty() {
                    egui_state.error_message = "Invalid SMILES".to_string();
                } else {
                    egui_state.error_message = e.clone();
                }
            }
        }
    }
}
