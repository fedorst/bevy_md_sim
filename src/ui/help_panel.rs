// src/ui/help_panel.rs

use super::UiSet;
use crate::components::Atom;
use crate::interaction::SelectionState;
use bevy::picking::hover::PickingInteraction;
use bevy::prelude::*;

#[derive(Component)]
struct HelpPanel;
#[derive(Resource, Default)]
struct HelpState {
    pub visible: bool,
}

pub struct HelpPanelPlugin;

impl Plugin for HelpPanelPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<HelpState>()
            .add_systems(Startup, setup_help_panel_ui)
            .add_systems(
                Update,
                (toggle_help_visibility, update_help_panel).in_set(UiSet),
            );
    }
}

fn setup_help_panel_ui(mut commands: Commands) {
    commands.spawn((
        Text::new("h: Toggle Help"),
        Node {
            position_type: PositionType::Absolute,
            bottom: Val::Px(10.0),
            right: Val::Px(10.0),
            padding: UiRect::all(Val::Px(10.0)),
            width: Val::Px(384.0),
            ..default()
        },
        BackgroundColor(Color::BLACK.with_alpha(0.75)),
        TextFont {
            font_size: 16.0,
            ..default()
        },
        HelpPanel,
    ));
}

// --- Systems (copied from original ui.rs) ---

fn toggle_help_visibility(keys: Res<ButtonInput<KeyCode>>, mut help_state: ResMut<HelpState>) {
    if keys.just_pressed(KeyCode::KeyH) {
        help_state.visible = !help_state.visible;
    }
}

fn update_help_panel(
    help_state: Res<HelpState>,
    selection: Res<SelectionState>,
    hover_query: Query<&PickingInteraction, With<Atom>>,
    mut help_panel_query: Query<(&mut Text, &mut Visibility), With<HelpPanel>>,
) {
    let Ok((mut text, mut visibility)) = help_panel_query.single_mut() else {
        return;
    };

    if !help_state.visible {
        *visibility = Visibility::Hidden;
        let prompt = "h: Toggle Help".to_string();
        if text.0 != prompt {
            text.0 = prompt;
        }
        return;
    }
    *visibility = Visibility::Visible;

    let is_hovering_atom = hover_query
        .iter()
        .any(|i| *i == PickingInteraction::Hovered);

    let mut lines = vec!["h: Toggle Help", "Space: Toggle Pause"];

    if is_hovering_atom {
        lines.push("LMB (click): Select / Deselect Atom");
    }

    lines.push("LMB (drag): Rotate Camera");
    lines.push("RMB (drag): Pan Camera");
    lines.push("F: Focus on last Selected");

    if !selection.selected.is_empty() {
        lines.push("Shift+LMB: Add/Remove from Selection");
        lines.push("Delete: Delete Selected");
    }

    if selection.selected.len() == 2 {
        lines.push("b: Toggle Bond");
    }

    text.0 = lines.join("\n");
}
