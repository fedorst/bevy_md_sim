// src/ui/pause_menu.rs

use super::UiSet;
use crate::resources::{
    ActiveInput, ForceMultiplier, PauseMenuState, SimulationParameters, SimulationState, Thermostat,
};
use bevy::input::keyboard::{Key, KeyboardInput};
use bevy::prelude::*;
use bevy::ui::widget::TextUiWriter;

#[derive(Component)]
struct PauseMenuPanel;
#[derive(Component)]
struct PauseMenuButton;
#[derive(Component)]
struct DtInputField;
#[derive(Component)]
struct TempInputField;
#[derive(Component)]
struct TauInputField;
#[derive(Component)]
struct ForceMultiplierInputField;
#[derive(Component)]
struct ApplyButton;

pub struct PauseMenuPlugin;

impl Plugin for PauseMenuPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_pause_menu_ui).add_systems(
            Update,
            (
                toggle_pause_menu,
                update_pause_menu_panel,
                handle_input_clicks,
                handle_keyboard_input,
                update_input_field_display,
                update_input_field_highlight,
                handle_apply_button,
            )
                .in_set(UiSet),
        );
    }
}

fn setup_pause_menu_ui(mut commands: Commands) {
    //  (m) Pause Menu Button
    commands
        .spawn((
            Button,
            PauseMenuButton,
            Node {
                position_type: PositionType::Absolute,
                top: Val::Px(10.0),
                right: Val::Px(10.0),
                padding: UiRect::all(Val::Px(10.0)),
                ..default()
            },
            BackgroundColor(Color::srgb(0.2, 0.2, 0.2)),
        ))
        .with_child((
            Text::new("(m) Pause Menu"),
            TextFont {
                font_size: 20.0,
                ..default()
            },
            TextColor(Color::WHITE),
        ));

    // --- The Pause Menu Panel (initially hidden) ---
    commands
        .spawn((
            PauseMenuPanel,
            Visibility::Hidden,
            Node {
                width: Val::Percent(100.0),
                height: Val::Percent(100.0),
                justify_content: JustifyContent::Center,
                align_items: AlignItems::Center,
                ..default()
            },
            BackgroundColor(Color::BLACK.with_alpha(0.75)), // Semi-transparent overlay
        ))
        .with_children(|parent| {
            // Main settings box
            parent
                .spawn(Node {
                    width: Val::Px(400.0),
                    flex_direction: FlexDirection::Column,
                    padding: UiRect::all(Val::Px(20.0)),
                    ..default()
                })
                .with_children(|p| {
                    // Title
                    p.spawn((
                        Text::new("Settings"),
                        TextFont {
                            font_size: 32.0,
                            ..default()
                        },
                        Node {
                            width: Val::Percent(100.0),
                            justify_content: JustifyContent::Center,
                            margin: UiRect::bottom(Val::Px(20.0)),
                            ..default()
                        },
                    ));

                    // Spawn the input rows
                    spawn_input_row(p, "Timestep (dt)", "0.0".to_string(), DtInputField);
                    spawn_input_row(p, "Target Temp (K)", "0.0".to_string(), TempInputField);
                    spawn_input_row(p, "Thermostat Tau (ps)", "0.0".to_string(), TauInputField);
                    spawn_input_row(
                        p,
                        "Force Multiplier",
                        "1.0".to_string(),
                        ForceMultiplierInputField,
                    );

                    // Apply Button
                    p.spawn((
                        Button,
                        ApplyButton,
                        Node {
                            margin: UiRect::top(Val::Px(20.0)),
                            padding: UiRect::all(Val::Px(10.0)),
                            ..default()
                        },
                        BackgroundColor(Color::srgb(0.2, 0.6, 0.2)),
                    ))
                    .with_child(Text::new("Apply & Resume"));
                });
        });
}

fn spawn_input_row(
    parent: &mut ChildSpawnerCommands,
    label: &str,
    initial_value: String,
    marker_component: impl Component,
) {
    parent
        .spawn(Node {
            width: Val::Percent(100.0),
            flex_direction: FlexDirection::Row,
            justify_content: JustifyContent::SpaceBetween,
            align_items: AlignItems::Center,
            margin: UiRect::bottom(Val::Px(10.0)),
            ..default()
        })
        .with_children(|p| {
            p.spawn(Text::new(label));
            p.spawn((
                Button,
                marker_component,
                Node {
                    width: Val::Px(150.0),
                    height: Val::Px(40.0),
                    justify_content: JustifyContent::Center,
                    align_items: AlignItems::Center,
                    border: UiRect::all(Val::Px(2.0)),
                    ..default()
                },
                BackgroundColor(Color::srgb(0.1, 0.1, 0.1)),
                BorderColor(Color::srgb(0.5, 0.5, 0.5)),
            ))
            .with_child(Text::new(initial_value));
        });
}

fn toggle_pause_menu(
    keys: Res<ButtonInput<KeyCode>>,
    mut sim_state: ResMut<SimulationState>,
    mut menu_state: ResMut<PauseMenuState>,
    interaction_query: Query<&Interaction, (Changed<Interaction>, With<PauseMenuButton>)>,
    sim_params: Res<SimulationParameters>,
    thermostat: Res<Thermostat>,
    force_multiplier: Res<ForceMultiplier>,
) {
    let mut should_toggle = keys.just_pressed(KeyCode::KeyM);
    if menu_state.visible && keys.just_pressed(KeyCode::Escape) {
        should_toggle = true;
    }
    if !should_toggle {
        for interaction in &interaction_query {
            if *interaction == Interaction::Pressed {
                should_toggle = true;
                break;
            }
        }
    }

    if should_toggle {
        if menu_state.visible {
            menu_state.visible = false;
        } else {
            menu_state.visible = true;
            sim_state.paused = true;
            menu_state.active_input = None;

            menu_state.dt_str = sim_params.dt.to_string();
            menu_state.temp_str = thermostat.target_temperature.to_string();
            menu_state.tau_str = thermostat.tau.to_string();
            menu_state.force_multiplier_str = force_multiplier.0.to_string();
        }
    }
}

fn update_pause_menu_panel(
    menu_state: Res<PauseMenuState>,
    mut panel_query: Query<&mut Visibility, With<PauseMenuPanel>>,
) {
    if menu_state.is_changed() {
        if let Ok(mut visibility) = panel_query.single_mut() {
            *visibility = if menu_state.visible {
                Visibility::Visible
            } else {
                Visibility::Hidden
            };
        }
    }
}

fn handle_input_clicks(
    mut menu_state: ResMut<PauseMenuState>,
    dt_q: Query<&Interaction, (Changed<Interaction>, With<DtInputField>)>,
    temp_q: Query<&Interaction, (Changed<Interaction>, With<TempInputField>)>,
    tau_q: Query<&Interaction, (Changed<Interaction>, With<TauInputField>)>,
    force_q: Query<&Interaction, (Changed<Interaction>, With<ForceMultiplierInputField>)>,
) {
    if let Ok(Interaction::Pressed) = dt_q.single() {
        menu_state.active_input = Some(ActiveInput::Dt);
    }
    if let Ok(Interaction::Pressed) = temp_q.single() {
        menu_state.active_input = Some(ActiveInput::Temp);
    }
    if let Ok(Interaction::Pressed) = tau_q.single() {
        menu_state.active_input = Some(ActiveInput::Tau);
    }
    if let Ok(Interaction::Pressed) = force_q.single() {
        menu_state.active_input = Some(ActiveInput::ForceMultiplier);
    }
}

fn handle_keyboard_input(
    mut menu_state: ResMut<PauseMenuState>,
    mut key_evr: EventReader<KeyboardInput>,
) {
    if let Some(active_input) = menu_state.active_input {
        let target_str = match active_input {
            ActiveInput::Dt => &mut menu_state.dt_str,
            ActiveInput::Temp => &mut menu_state.temp_str,
            ActiveInput::Tau => &mut menu_state.tau_str,
            ActiveInput::ForceMultiplier => &mut menu_state.force_multiplier_str,
        };

        for ev in key_evr.read() {
            if !ev.state.is_pressed() {
                continue;
            }
            match &ev.logical_key {
                Key::Character(chars) => {
                    for char in chars.chars() {
                        if char.is_ascii_digit() || (char == '.' && !target_str.contains('.')) {
                            target_str.push(char);
                        }
                    }
                }
                Key::Backspace => {
                    target_str.pop();
                }
                _ => {}
            }
        }
    }
}

fn update_input_field_display(
    menu_state: Res<PauseMenuState>,
    mut writer: TextUiWriter,
    dt_q: Query<&Children, With<DtInputField>>,
    temp_q: Query<&Children, With<TempInputField>>,
    tau_q: Query<&Children, With<TauInputField>>,
    force_q: Query<&Children, With<ForceMultiplierInputField>>,
) {
    if menu_state.is_changed() {
        if let Ok(children) = dt_q.single() {
            if let Some(text_entity) = children.first() {
                *writer.text(*text_entity, 0) = menu_state.dt_str.clone();
            }
        }
        if let Ok(children) = temp_q.single() {
            if let Some(text_entity) = children.first() {
                *writer.text(*text_entity, 0) = menu_state.temp_str.clone();
            }
        }
        if let Ok(children) = tau_q.single() {
            if let Some(text_entity) = children.first() {
                *writer.text(*text_entity, 0) = menu_state.tau_str.clone();
            }
        }
        if let Ok(children) = force_q.single() {
            if let Some(text_entity) = children.first() {
                *writer.text(*text_entity, 0) = menu_state.force_multiplier_str.clone();
            }
        }
    }
}

fn update_input_field_highlight(
    mut commands: Commands,
    menu_state: Res<PauseMenuState>,
    dt_e: Query<Entity, With<DtInputField>>,
    temp_e: Query<Entity, With<TempInputField>>,
    tau_e: Query<Entity, With<TauInputField>>,
    force_e: Query<Entity, With<ForceMultiplierInputField>>,
) {
    let mut clear_outline = |entity: Entity| {
        if let Ok(mut entity_commands) = commands.get_entity(entity) {
            entity_commands.remove::<Outline>();
        }
    };
    if let Ok(entity) = dt_e.single() {
        clear_outline(entity);
    }
    if let Ok(entity) = temp_e.single() {
        clear_outline(entity);
    }
    if let Ok(entity) = tau_e.single() {
        clear_outline(entity);
    }
    if let Ok(entity) = force_e.single() {
        clear_outline(entity);
    }
    if let Some(active_input) = menu_state.active_input {
        let target_entity_result = match active_input {
            ActiveInput::Dt => dt_e.single(),
            ActiveInput::Temp => temp_e.single(),
            ActiveInput::Tau => tau_e.single(),
            ActiveInput::ForceMultiplier => force_e.single(),
        };
        if let Ok(target_entity) = target_entity_result {
            commands.entity(target_entity).insert(Outline {
                width: Val::Px(2.0),
                offset: Val::Px(2.0),
                color: Color::WHITE,
            });
        }
    }
}

fn handle_apply_button(
    interaction_q: Query<&Interaction, (Changed<Interaction>, With<ApplyButton>)>,
    mut menu_state: ResMut<PauseMenuState>,
    mut sim_state: ResMut<SimulationState>,
    mut sim_params: ResMut<SimulationParameters>,
    mut thermostat: ResMut<Thermostat>,
    mut force_multiplier: ResMut<ForceMultiplier>,
) {
    if let Ok(Interaction::Pressed) = interaction_q.single() {
        if let Ok(dt) = menu_state.dt_str.parse::<f32>() {
            sim_params.dt = dt;
        }
        if let Ok(temp) = menu_state.temp_str.parse::<f32>() {
            thermostat.target_temperature = temp;
        }
        if let Ok(tau) = menu_state.tau_str.parse::<f32>() {
            thermostat.tau = tau;
        }
        if let Ok(fm) = menu_state.force_multiplier_str.parse::<f32>() {
            force_multiplier.0 = fm;
        }

        menu_state.visible = false;
        sim_state.paused = false;
    }
}
