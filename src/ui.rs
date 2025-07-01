use super::resources::*;
use super::setup::{EnergyDisplayText, PauseText, TempDisplayText, TimeDisplayText}; // Use the public marker component
use bevy::prelude::*;

pub struct UIPlugin;

impl Plugin for UIPlugin {
    fn build(&self, app: &mut App) {
        // We can run these systems unconditionally
        app.add_systems(
            Update,
            (
                handle_simulation_control,
                continuous_simulation,
                update_pause_text,
                update_time_display,
                track_active_wall_time,
                update_temp_display,
                update_energy_display,
            ),
        );
    }
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

fn track_active_wall_time(
    time: Res<Time>,
    sim_state: Res<SimulationState>,
    mut active_wall_time: ResMut<ActiveWallTime>,
) {
    if !sim_state.paused {
        active_wall_time.0 += time.delta_secs();
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

fn handle_simulation_control(
    mut state: ResMut<SimulationState>,
    mut step: ResMut<StepSimulation>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    if keys.just_pressed(KeyCode::Space) {
        state.paused = !state.paused;
    }
    // Right arrow steps one frame when paused
    if state.paused && keys.just_pressed(KeyCode::ArrowRight) {
        step.0 = true;
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

fn continuous_simulation(mut step: ResMut<StepSimulation>, state: Res<SimulationState>) {
    if !state.paused {
        step.0 = true;
    }
}
