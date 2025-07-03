// src/ui/mod.rs

mod help_panel;
mod hud;
mod info_panel;
mod pause_menu;

use crate::interaction::InteractionSet;
use crate::resources::{ActiveWallTime, PauseMenuState, SimulationState, StepSimulation};
use bevy::prelude::*;

use help_panel::HelpPanelPlugin;
use hud::HudPlugin;
use info_panel::InfoPanelPlugin;
use pause_menu::PauseMenuPlugin;

// This set can be used if any UI systems need to be ordered relative to each other.
#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct UiSet;

pub struct UIPlugin;

impl Plugin for UIPlugin {
    fn build(&self, app: &mut App) {
        app
            // Add any resources or events shared across all UI modules here
            .configure_sets(Update, UiSet.after(InteractionSet))
            .add_plugins((HudPlugin, InfoPanelPlugin, HelpPanelPlugin, PauseMenuPlugin))
            // Keep global UI controls here
            .add_systems(
                Update,
                (
                    handle_simulation_control,
                    continuous_simulation,
                    track_active_wall_time,
                )
                    .in_set(UiSet),
            );
    }
}

// These systems are global controls, not tied to a specific panel,
// so they can live in the main `mod.rs`.

fn handle_simulation_control(
    mut sim_state: ResMut<SimulationState>,
    mut menu_state: ResMut<PauseMenuState>,
    mut step: ResMut<StepSimulation>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    if keys.just_pressed(KeyCode::Space) {
        if menu_state.visible {
            // If the menu is open, Space closes it and ALWAYS unpauses.
            menu_state.visible = false;
            sim_state.paused = false;
        } else {
            // Otherwise, it's just a simple toggle.
            sim_state.paused = !sim_state.paused;
        }
    }
    // Right arrow steps one frame when paused
    if sim_state.paused && keys.just_pressed(KeyCode::ArrowRight) {
        step.0 = true;
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

fn continuous_simulation(mut step: ResMut<StepSimulation>, state: Res<SimulationState>) {
    if !state.paused {
        step.0 = true;
    }
}
