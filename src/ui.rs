use super::resources::*;
use crate::components::{Atom, Force, Velocity};
use crate::interaction::InteractionSet; // Import the InteractionSet
use crate::interaction::SelectionState; // Import from the new location
use crate::resources::BondOrder;
use bevy::picking::hover::PickingInteraction;
use bevy::prelude::*;
use bevy::ui::widget::TextUiWriter;

#[derive(Component)]
pub struct PauseText;

#[derive(Component)]
pub struct EnergyDisplayText;

#[derive(Component)]
pub struct TempDisplayText;

#[derive(Component)]
pub struct TimeDisplayText;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct UiSet;

pub struct UIPlugin;

#[derive(Component)]
pub struct HelpPanel;

#[derive(Resource, Default)]
pub struct HelpState {
    pub visible: bool,
}

#[derive(Component)]
struct PauseMenuPanel;

#[derive(Component)]
struct PauseMenuButton;

#[derive(Component)]
pub struct DebugInfoPanel;

impl Plugin for UIPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<HelpState>()
            .init_resource::<PauseMenuState>()
            .add_systems(Startup, setup_ui)
            .add_systems(
                Update,
                (
                    handle_simulation_control,
                    toggle_pause_menu,
                    continuous_simulation,
                    update_pause_text,
                    update_time_display,
                    track_active_wall_time,
                    update_temp_display,
                    update_energy_display,
                    update_info_panel,
                    update_pause_menu_panel,
                    toggle_help_visibility,
                    update_help_panel,
                )
                    .in_set(UiSet)
                    .after(InteractionSet),
            );
    }
}

fn toggle_pause_menu(
    keys: Res<ButtonInput<KeyCode>>,
    mut sim_state: ResMut<SimulationState>,
    mut menu_state: ResMut<PauseMenuState>,
    interaction_query: Query<&Interaction, (Changed<Interaction>, With<PauseMenuButton>)>,
) {
    let mut should_toggle = keys.just_pressed(KeyCode::KeyM);
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
            // If menu is open, just close it. Do not touch sim_state.
            menu_state.visible = false;
        } else {
            // If menu is closed, open it AND pause the simulation.
            menu_state.visible = true;
            sim_state.paused = true;
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
        // If help is hidden, just show the prompt.
        text.0 = "h: Toggle Help".to_string();
        return;
    }

    // 1. Control overall visibility first
    if !help_state.visible {
        *visibility = Visibility::Hidden;
        if !text.0.is_empty() {
            text.0.clear();
        }
        return; // Don't do any work if it's hidden
    }
    *visibility = Visibility::Visible;

    // 2. Check for hover state
    let is_hovering_atom = hover_query
        .iter()
        .any(|i| *i == PickingInteraction::Hovered);

    // 3. Build the text string piece by piece
    let mut lines = vec!["h: Toggle Help", "Space: Toggle Pause"];

    // 4. Add context-specific lines
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

fn display_single_atom_info(
    entity: Entity,
    atom_id_map: &Res<AtomIdMap>,
    atom_query: &Query<(&Atom, &Transform, &Velocity, &Force)>,
) -> String {
    if let Ok((atom, _transform, velocity, force)) = atom_query.get(entity) {
        let pretty_id = atom_id_map.entity_to_id.get(&entity).unwrap();
        format!(
            "Selected: {}\n\
             Type: {}\n\
             Vel (mag): {:.2}\n\
             --- Forces (mag) ---\n\
             Total: {:.1}\n\
             Bond: {:.1}\n\
             Angle: {:.1}\n\
             Non-Bonded: {:.1}",
            pretty_id,
            atom.type_name,
            velocity.length(),
            force.total_magnitude(),
            force.bond.length(),
            force.angle.length(),
            force.non_bonded.length()
        )
    } else {
        "Error: Could not get atom data.".to_string()
    }
}

fn display_bond_info(
    e1: Entity,
    e2: Entity,
    connectivity: &Res<SystemConnectivity>,
    force_field: &Res<ForceField>,
    atom_id_map: &Res<AtomIdMap>,
    atom_query: &Query<(&Atom, &Transform, &Velocity, &Force)>,
) -> Option<String> {
    // Return Option<String> in case there's no bond
    if let Some(bond) = connectivity
        .bonds
        .iter()
        .find(|b| (b.a == e1 && b.b == e2) || (b.a == e2 && b.b == e1))
    {
        if let Ok([(type1, t1, _, _), (type2, t2, _, _)]) = atom_query.get_many([e1, e2]) {
            if let Some(&(_k, r0)) = force_field.bond_params.get(&(
                type1.type_name.clone(),
                type2.type_name.clone(),
                bond.order,
            )) {
                let current_dist = t1.translation.distance(t2.translation);
                let id1 = atom_id_map.entity_to_id.get(&e1).unwrap();
                let id2 = atom_id_map.entity_to_id.get(&e2).unwrap();
                let bond_type_str = match bond.order {
                    BondOrder::Single => "Single Bond",
                    BondOrder::Double => "Double Bond",
                };
                return Some(format!(
                    "Selected {}-{}\n\
                    Type: {} ({}-{})\n\
                     --- Geometry ---\n\
                     Current Len: {:.3} nm\n\
                     Optimal Len: {:.3} nm\n\
                     Strain: {:.2}%",
                    id1,
                    id2,
                    bond_type_str,
                    type1.type_name,
                    type2.type_name,
                    current_dist,
                    r0,
                    (current_dist / r0 - 1.0) * 100.0
                ));
            }
        }
    }
    None // Return None if no bond was found
}

fn display_angle_info(
    entities: &[Entity],
    connectivity: &Res<SystemConnectivity>,
    force_field: &Res<ForceField>,
    atom_id_map: &Res<AtomIdMap>,
    atom_query: &Query<(&Atom, &Transform, &Velocity, &Force)>,
) -> Option<String> {
    if let Some(angle) = connectivity.angles.iter().find(|a| {
        entities.contains(&a.a) && entities.contains(&a.center) && entities.contains(&a.b)
    }) {
        if let Ok(
            [
                (type_a, t_a, _, _),
                (type_center, t_center, _, _),
                (type_b, t_b, _, _),
            ],
        ) = atom_query.get_many([angle.a, angle.center, angle.b])
        {
            if let Some(&(_k, theta0)) = force_field.angle_params.get(&(
                type_a.type_name.clone(),
                type_center.type_name.clone(),
                type_b.type_name.clone(),
            )) {
                let v1 = t_a.translation - t_center.translation;
                let v2 = t_b.translation - t_center.translation;
                let current_angle_rad = v1.angle_between(v2);
                let id_a = atom_id_map.entity_to_id.get(&angle.a).unwrap();
                let id_center = atom_id_map.entity_to_id.get(&angle.center).unwrap();
                let id_b = atom_id_map.entity_to_id.get(&angle.b).unwrap();
                return Some(format!(
                    "Selected Angle: {}-{}-{}\n\
                     --- Geometry ---\n\
                     Optimal: {:.2}deg\n\
                     Current: {:.2}deg ({:+.2})",
                    id_a,
                    id_center,
                    id_b,
                    theta0.to_degrees(),
                    current_angle_rad.to_degrees(),
                    current_angle_rad.to_degrees() - theta0.to_degrees()
                ));
            }
        }
    }
    None
}

fn update_info_panel(
    selection: Res<SelectionState>,
    connectivity: Res<SystemConnectivity>,
    force_field: Res<ForceField>,
    atom_id_map: Res<AtomIdMap>,
    atom_query: Query<(&Atom, &Transform, &Velocity, &Force)>,
    panel_query: Query<Entity, With<DebugInfoPanel>>,
    mut writer: TextUiWriter,
) {
    if let Ok(panel_entity) = panel_query.single() {
        let info_text = match selection.selected.len() {
            0 => "Select an atom for details".to_string(),
            1 => display_single_atom_info(selection.selected[0], &atom_id_map, &atom_query),
            2 => display_bond_info(
                selection.selected[0],
                selection.selected[1],
                &connectivity,
                &force_field,
                &atom_id_map,
                &atom_query,
            )
            .unwrap_or_else(|| format!("Selected {} atoms.", selection.selected.len())),
            3 => display_angle_info(
                &selection.selected,
                &connectivity,
                &force_field,
                &atom_id_map,
                &atom_query,
            )
            .unwrap_or_else(|| format!("Selected {} atoms.", selection.selected.len())),
            _ => format!("Selected {} atoms.", selection.selected.len()),
        };
        *writer.text(panel_entity, 0) = info_text;
    }
}

fn setup_ui(mut commands: Commands) {
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
    commands.spawn((
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
    ));

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
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            bottom: Val::Px(10.0),
            left: Val::Px(10.0),
            width: Val::Px(384.0),
            padding: UiRect::all(Val::Px(10.0)),
            flex_direction: FlexDirection::Column,
            ..default()
        },
        BackgroundColor(Color::BLACK.with_alpha(0.75)),
        Text::new("Select an atom for details"),
        TextFont {
            font_size: 16.0,
            ..default()
        },
        DebugInfoPanel,
    ));
    commands.spawn((
        // It starts with the default "collapsed" text
        Text::new("h: Toggle Help"),
        Node {
            position_type: PositionType::Absolute,
            bottom: Val::Px(10.0),
            right: Val::Px(10.0), // Let's put it in the bottom-right
            padding: UiRect::all(Val::Px(10.0)),
            width: Val::Px(384.0),
            ..default()
        },
        BackgroundColor(Color::BLACK.with_alpha(0.75)),
        TextFont {
            font_size: 16.0,
            ..default()
        },
        // Add our marker component so we can find it.
        HelpPanel,
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
