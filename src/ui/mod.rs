// src/ui/mod.rs

mod force_inspector;
mod help_panel;
mod hud;
mod info_panel;
mod pause_menu;

use crate::interaction::InteractionSet;
use crate::resources::{SimulationState, StepSimulation};
use bevy::prelude::*;

use force_inspector::ForceInspectorPlugin;
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
        app.configure_sets(Update, UiSet.after(InteractionSet))
            .add_plugins((
                HudPlugin,
                InfoPanelPlugin,
                HelpPanelPlugin,
                PauseMenuPlugin,
                ForceInspectorPlugin,
            ))
            // Keep global UI controls here
            .add_systems(
                Update,
                (handle_simulation_control, continuous_simulation).in_set(UiSet),
            );
    }
}

// These systems are global controls, not tied to a specific panel,
// so they can live in the main `mod.rs`.

fn handle_simulation_control(
    mut sim_state: ResMut<SimulationState>,
    mut step: ResMut<StepSimulation>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    // Space bar ONLY toggles the paused state of the simulation.
    if keys.just_pressed(KeyCode::Space) {
        sim_state.paused = !sim_state.paused;
    }

    // Right arrow still steps one frame when paused.
    if sim_state.paused && keys.just_pressed(KeyCode::ArrowRight) {
        step.0 = true;
    }
}

fn continuous_simulation(mut step: ResMut<StepSimulation>, state: Res<SimulationState>) {
    if !state.paused {
        step.0 = true;
    }
}
