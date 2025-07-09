// src/ui/hud.rs

use crate::resources::{
    ActiveWallTime, CurrentTemperature, SimulationParameters, SimulationState, StepCount,
    SystemEnergy, Thermostat, ThermostatScale,
};
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

pub struct HudPlugin;

impl Plugin for HudPlugin {
    fn build(&self, app: &mut App) {
        // The plugin now only needs to register one system.
        app.add_systems(EguiPrimaryContextPass, hud_egui_system);
    }
}

/// This single system draws all the heads-up display elements using egui.
fn hud_egui_system(
    mut contexts: EguiContexts,
    // Query for all the data resources needed for the display.
    sim_state: Res<SimulationState>,
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

                // Use a rich text format for a nice header
                ui.label(
                    egui::RichText::new("Simulation Info").font(egui::FontId::proportional(16.0)),
                );
                ui.separator();

                // Build the text content dynamically
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
    if sim_state.paused {
        egui::Area::new(egui::Id::new("hud_paused_text"))
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, -100.0)) // Position slightly above center
            .show(ctx, |ui| {
                // Use a large, bold, semi-transparent font for the status.
                let text = egui::RichText::new("PAUSED")
                    .font(egui::FontId::proportional(48.0))
                    .color(egui::Color32::from_white_alpha(180))
                    .strong();
                ui.label(text);
            });
    }
}
