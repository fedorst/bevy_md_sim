// src/ui/hud.rs

use super::UiSet;
use crate::resources::{
    ActiveWallTime, CurrentTemperature, SimulationParameters, SimulationState, StepCount,
    SystemEnergy, Thermostat, ThermostatScale,
};
use crate::spawning::{SMILESValidationResult, SpawnMoleculeFromSMILESEvent, ValidateSMILESEvent};
use bevy::prelude::*;
use bevy_egui::EguiPrimaryContextPass;
use bevy_egui::{EguiContexts, egui};

#[derive(Resource, Default)]
struct EguiInputState {
    smiles: String,
    is_valid: bool,
    error_message: String,
}

#[derive(Component)]
struct EnergyDisplayText;
#[derive(Component)]
struct TempDisplayText;
#[derive(Component)]
struct TimeDisplayText;
#[derive(Component)]
struct PauseText;

pub struct HudPlugin;

impl Plugin for HudPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<EguiInputState>()
            .add_systems(Startup, setup_hud_ui)
            .add_systems(EguiPrimaryContextPass, smiles_input_egui_system)
            .add_systems(
                Update,
                (
                    update_validation_state,
                    update_pause_text,
                    update_time_display,
                    update_energy_display,
                    update_temp_display,
                )
                    .after(UiSet),
            );
    }
}
fn smiles_input_egui_system(
    mut contexts: EguiContexts,
    mut egui_state: ResMut<EguiInputState>,
    mut validate_writer: EventWriter<ValidateSMILESEvent>,
    mut spawn_writer: EventWriter<SpawnMoleculeFromSMILESEvent>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    let window = egui::Window::new("Molecule Input")
        .anchor(egui::Align2::CENTER_BOTTOM, egui::vec2(0.0, -20.0))
        .collapsible(false)
        .resizable(false)
        .title_bar(false);

    // --- THE FIX: Use the resolved `ctx` to get the style ---
    let frame_stroke = if egui_state.is_valid {
        egui::Stroke::new(1.0, egui::Color32::from_rgb(0, 255, 255)) // Cyan
    } else {
        egui::Stroke::new(2.0, egui::Color32::from_rgb(255, 50, 50)) // Red
    };
    let frame = egui::Frame::group(&ctx.style()).stroke(frame_stroke);

    // --- THE FIX: Pass the resolved `ctx` to the show method ---
    window.frame(frame).show(ctx, |ui| {
        ui.set_width(300.0);
        ui.label("SMILES String:");

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
                spawn_writer.write(SpawnMoleculeFromSMILESEvent(egui_state.smiles.clone()));
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

fn setup_hud_ui(mut commands: Commands) {
    commands.spawn((
        Text::new("Time: ..."),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(10.0),
            left: Val::Px(10.0),
            ..default()
        },
        TextFont {
            font_size: 20.0,
            ..default()
        },
        TextColor(Color::WHITE),
        TimeDisplayText,
    ));
    commands.spawn((
        Text::new("PAUSED"),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(10.0),
            right: Val::Px(10.0),
            left: Val::Percent(50.0),
            ..default()
        },
        Transform::from_translation(Vec3::new(-50.0, 0.0, 0.0)),
        TextFont {
            font_size: 30.0,
            ..default()
        },
        TextColor(Color::WHITE),
        PauseText,
    ));

    commands.spawn((
        Text::new("Temp: ..."),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(80.0),
            left: Val::Px(10.0),
            ..default()
        },
        TextFont {
            font_size: 20.0,
            ..default()
        },
        TextColor(Color::WHITE),
        TempDisplayText,
    ));

    commands.spawn((
        Text::new("Energy: ..."),
        Node {
            position_type: PositionType::Absolute,
            // Position it below the temperature display
            top: Val::Px(130.0),
            left: Val::Px(10.0),
            width: Val::Px(384.0),
            ..default()
        },
        TextFont {
            font_size: 16.0,
            ..default()
        },
        TextColor(Color::WHITE),
        EnergyDisplayText,
    ));
}

fn update_energy_display(
    energy: Res<SystemEnergy>,
    mut query: Query<&mut Text, With<EnergyDisplayText>>,
) {
    if let Ok(mut text) = query.single_mut() {
        text.0 = format!(
            "--- Energy (kJ/mol) ---\n\
             Total: {:.2}\n\
             Potential: {:.2}\n\
             Kinetic: {:.2}",
            energy.total, energy.potential, energy.kinetic
        );
    }
}

fn update_temp_display(
    thermostat: Res<Thermostat>,
    current_temp: Res<CurrentTemperature>,
    thermostat_scale: Res<ThermostatScale>,
    mut query: Query<&mut Text, With<TempDisplayText>>,
) {
    if let Ok(mut text) = query.single_mut() {
        text.0 = format!(
            "Temp: {:.1} K / {:.1} K\nScale: {:.6}",
            current_temp.0, thermostat.target_temperature, thermostat_scale.0
        );
    }
}

fn update_time_display(
    active_wall_time: Res<ActiveWallTime>,
    step_count: Res<StepCount>,
    params: Res<SimulationParameters>,
    mut query: Query<&mut Text, With<TimeDisplayText>>,
) {
    if let Ok(mut text) = query.single_mut() {
        // `dt` is in picoseconds (ps)?
        let simulated_time_ps = step_count.0 as f32 * params.dt;
        text.0 = format!(
            "Sim Time: {:.3} ps\nWall Time: {:.2} s\nSteps: {:}",
            simulated_time_ps, active_wall_time.0, step_count.0
        );
    }
}

fn update_pause_text(state: Res<SimulationState>, mut query: Query<&mut Text, With<PauseText>>) {
    if state.is_changed() {
        if let Ok(mut text) = query.single_mut() {
            text.0 = if state.paused {
                "PAUSED".to_string()
            } else {
                "RUNNING".to_string()
            };
        }
    }
}
