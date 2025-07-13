// src/ui/mod.rs

mod force_inspector;
mod help_panel;
mod hud;
mod info_panel;
mod pause_menu;
pub mod peptide_builder;
mod system_metrics_panel;

use crate::AppState;
use crate::interaction::InteractionSet;
use bevy::prelude::*;
use force_inspector::ForceInspectorPlugin;
use help_panel::HelpPanelPlugin;
use hud::HudPlugin;
use info_panel::InfoPanelPlugin;
use pause_menu::PauseMenuPlugin;
use peptide_builder::PeptideBuilderPlugin;
use system_metrics_panel::SystemMetricsPanelPlugin;
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
                SystemMetricsPanelPlugin,
                PeptideBuilderPlugin,
            ))
            // Keep global UI controls here
            .add_systems(Update, (handle_simulation_control).in_set(UiSet));
    }
}

fn handle_simulation_control(
    app_state: Res<State<AppState>>,
    mut next_state: ResMut<NextState<AppState>>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    // Space bar toggles between Running and Paused.
    if keys.just_pressed(KeyCode::Space) {
        match app_state.get() {
            AppState::Running => next_state.set(AppState::Paused),
            AppState::Paused => next_state.set(AppState::Running),
            _ => {} // Do nothing while initializing
        }
    }
}
