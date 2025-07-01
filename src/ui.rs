use super::resources::*;
use crate::components::{AtomType, Force, Velocity};
use crate::interaction::InteractionSet; // Import the InteractionSet
use crate::interaction::SelectionState; // Import from the new location

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
pub struct DebugInfoPanel;

impl Plugin for UIPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_ui).add_systems(
            Update,
            (
                handle_simulation_control,
                continuous_simulation,
                update_pause_text,
                update_time_display,
                track_active_wall_time,
                update_temp_display,
                update_energy_display,
                update_info_panel,
            )
                .in_set(UiSet)
                .after(InteractionSet),
        );
    }
}

fn display_single_atom_info(
    entity: Entity,
    atom_id_map: &Res<AtomIdMap>,
    atom_query: &Query<(&AtomType, &Transform, &Velocity, &Force)>,
) -> String {
    if let Ok((atom_type, _transform, velocity, force)) = atom_query.get(entity) {
        let pretty_id = atom_id_map.entity_to_id.get(&entity).unwrap();
        format!(
            "Selected: {}\n\
             Type: {:?}\n\
             Vel (mag): {:.2}\n\
             --- Forces (mag) ---\n\
             Total: {:.1}\n\
             Bond: {:.1}\n\
             Angle: {:.1}\n\
             Non-Bonded: {:.1}",
            pretty_id,
            atom_type,
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
    atom_query: &Query<(&AtomType, &Transform, &Velocity, &Force)>,
) -> Option<String> {
    // Return Option<String> in case there's no bond
    if let Some(_bond) = connectivity
        .bonds
        .iter()
        .find(|b| (b.a == e1 && b.b == e2) || (b.a == e2 && b.b == e1))
    {
        if let Ok([(type1, t1, _, _), (type2, t2, _, _)]) = atom_query.get_many([e1, e2]) {
            if let Some(&(_k, r0)) = force_field.bond_params.get(&(*type1, *type2)) {
                let current_dist = t1.translation.distance(t2.translation);
                let id1 = atom_id_map.entity_to_id.get(&e1).unwrap();
                let id2 = atom_id_map.entity_to_id.get(&e2).unwrap();
                return Some(format!(
                    "Selected Bond:\n{}-{}\n\
                     --- Geometry ---\n\
                     Current Len: {:.3} nm\n\
                     Optimal Len: {:.3} nm\n\
                     Strain: {:.2}%",
                    id1,
                    id2,
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
    atom_query: &Query<(&AtomType, &Transform, &Velocity, &Force)>,
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
            if let Some(&(_k, theta0)) =
                force_field
                    .angle_params
                    .get(&(*type_a, *type_center, *type_b))
            {
                let v1 = t_a.translation - t_center.translation;
                let v2 = t_b.translation - t_center.translation;
                let current_angle_rad = v1.angle_between(v2);
                let id_a = atom_id_map.entity_to_id.get(&angle.a).unwrap();
                let id_center = atom_id_map.entity_to_id.get(&angle.center).unwrap();
                let id_b = atom_id_map.entity_to_id.get(&angle.b).unwrap();
                return Some(format!(
                    "Selected Angle:\n{}-{}-{}\n\
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
    atom_query: Query<(&AtomType, &Transform, &Velocity, &Force)>,
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
            ..default()
        },
        TextFont {
            font_size: 30.0,
            ..default()
        },
        TextColor(Color::WHITE),
        PauseText,
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
            width: Val::Px(250.0),
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
            width: Val::Px(250.0),
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
