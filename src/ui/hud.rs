// src/ui/hud.rs

use super::UiSet;
use crate::resources::{
    ActiveWallTime, CurrentTemperature, SimulationParameters, SimulationState, StepCount,
    SystemEnergy, Thermostat, ThermostatScale,
};
use bevy::prelude::*;

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
        app.add_systems(Startup, setup_hud_ui).add_systems(
            Update,
            (
                update_pause_text,
                update_time_display,
                update_temp_display,
                update_energy_display,
            )
                .in_set(UiSet),
        );
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
