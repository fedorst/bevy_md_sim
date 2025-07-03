// src/ui/info_panel.rs

use super::UiSet;
use crate::components::{Atom, Force, Velocity};
use crate::interaction::SelectionState;
use crate::resources::{AtomIdMap, BondOrder, ForceField, SystemConnectivity};
use bevy::prelude::*;
use bevy::ui::widget::TextUiWriter;

#[derive(Component)]
pub struct DebugInfoPanel;

pub struct InfoPanelPlugin;

impl Plugin for InfoPanelPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_info_panel_ui)
            .add_systems(Update, update_info_panel.in_set(UiSet));
    }
}

fn setup_info_panel_ui(mut commands: Commands) {
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
}

// --- Systems and Helpers (copied from original ui.rs) ---

fn update_info_panel(
    selection: Res<SelectionState>,
    connectivity: Res<SystemConnectivity>,
    force_field: Res<ForceField>,
    atom_id_map: Res<AtomIdMap>,
    atom_query: Query<(&Atom, &Transform, &Velocity, &Force)>,
    panel_query: Query<Entity, With<DebugInfoPanel>>,
    mut writer: TextUiWriter,
) {
    let Ok(panel_entity) = panel_query.single() else {
        return;
    };
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
    None
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
