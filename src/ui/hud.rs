// src/ui/hud.rs

use super::UiSet;
use crate::resources::{
    ActiveWallTime, CurrentTemperature, SimulationParameters, SimulationState, StepCount,
    SystemEnergy, Thermostat, ThermostatScale,
};
use crate::spawning::{SMILESValidationResult, SpawnMoleculeFromSMILESEvent, ValidateSMILESEvent};
use bevy::input::keyboard::{Key, KeyboardInput};
use bevy::prelude::*;
use bevy::ui::widget::TextUiWriter;

// Add this new resource
#[derive(Resource, Default)]
struct SmilesInput {
    text: String,
    active: bool,
}

// Add this new component
#[derive(Component)]
struct SmilesInputField;

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
        app.init_resource::<SmilesInput>()
            .add_systems(Startup, (setup_hud_ui, setup_smiles_input_ui))
            .add_systems(
                Update,
                (
                    update_pause_text,
                    update_time_display,
                    update_temp_display,
                    update_energy_display,
                    handle_smiles_input_click,
                    handle_smiles_keyboard,
                    update_smiles_display,
                    handle_validation_ui,
                )
                    .in_set(UiSet),
            );
    }
}

fn setup_smiles_input_ui(mut commands: Commands) {
    commands
        .spawn((
            Button,
            SmilesInputField,
            Node {
                position_type: PositionType::Absolute,
                bottom: Val::Px(10.0),
                left: Val::Percent(50.0),
                width: Val::Px(300.0),
                border: UiRect::all(Val::Px(2.0)),
                padding: UiRect::all(Val::Px(5.0)),
                justify_content: JustifyContent::FlexStart,
                align_items: AlignItems::Center,
                ..default()
            },
            Transform::from_translation(Vec3::new(-150.0, 0.0, 0.0)),
            BackgroundColor(Color::srgb(0.1, 0.1, 0.1)),
            BorderColor(Color::WHITE),
        ))
        .with_child(Text::new("Click to type SMILES..."));
}

fn handle_smiles_input_click(
    mut interaction_q: Query<
        (&Interaction, &mut BorderColor),
        (Changed<Interaction>, With<SmilesInputField>),
    >,
    mut input_state: ResMut<SmilesInput>,
) {
    if let Ok((interaction, mut border)) = interaction_q.single_mut() {
        if *interaction == Interaction::Pressed {
            input_state.active = true;
            *border = Color::WHITE.into();
        }
    }
}

fn update_smiles_display(
    state: Res<SmilesInput>,
    input_field_q: Query<&Children, With<SmilesInputField>>,
    mut text_writer: TextUiWriter,
) {
    if state.is_changed() {
        if let Ok(children) = input_field_q.single() {
            if let Some(text_entity) = children.first() {
                *text_writer.text(*text_entity, 0) = if state.text.is_empty() && !state.active {
                    "Click to type SMILES...".to_string()
                } else {
                    format!("{}{}", state.text, if state.active { "_" } else { "" })
                };
            }
        }
    }
}

fn handle_validation_ui(
    mut validation_results: EventReader<SMILESValidationResult>,
    mut border_q: Query<&mut BorderColor, With<SmilesInputField>>,
    input_state: Res<SmilesInput>,
) {
    if let Ok(mut border) = border_q.single_mut() {
        if let Some(event) = validation_results.read().last() {
            info!("[UI] Received validation result: {:?}", event.0);
            if input_state.active {
                match &event.0 {
                    Ok(_) => *border = Color::WHITE.into(),
                    Err(e) if !e.trim().is_empty() => {
                        *border = Color::linear_rgba(0.8, 0.1, 0.1, 1.0).into()
                    }
                    Err(_) => *border = Color::linear_rgba(0.8, 0.1, 0.1, 1.0).into(),
                }
            }
        }

        if input_state.is_changed() && !input_state.active {
            *border = Color::WHITE.into();
        }
    }
}

fn handle_smiles_keyboard(
    mut input_state: ResMut<SmilesInput>,
    mut key_evr: EventReader<KeyboardInput>,
    // --- CHANGE these writers ---
    mut spawn_writer: EventWriter<SpawnMoleculeFromSMILESEvent>,
    mut validate_writer: EventWriter<ValidateSMILESEvent>,
) {
    if !input_state.active {
        return;
    }

    let mut string_changed = false;
    for ev in key_evr.read() {
        if !ev.state.is_pressed() {
            continue;
        }

        match &ev.logical_key {
            Key::Character(chars) => {
                input_state.text.push_str(chars);
                string_changed = true;
            }
            Key::Backspace => {
                input_state.text.pop();
                string_changed = true;
            }
            Key::Enter => {
                if !input_state.text.is_empty() {
                    // On Enter, send the final spawn event
                    spawn_writer.write(SpawnMoleculeFromSMILESEvent(input_state.text.clone()));
                }
                input_state.text.clear();
                input_state.active = false;
                string_changed = true; // To reset color
            }
            Key::Escape => {
                input_state.active = false;
                string_changed = true; // To reset color
            }
            _ => {}
        }
    }

    // If the string changed, send a validation event.
    if string_changed && input_state.active {
        info!(
            "[UI] String changed. Sending validation event for: '{}'",
            input_state.text
        );
        validate_writer.write(ValidateSMILESEvent(input_state.text.clone()));
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

// Move these systems from the old ui.rs into here.
// Their code does not need to change, but you'll need to update their `use` statements.

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
