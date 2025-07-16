use crate::AppState;
use crate::components::{Acceleration, Atom, Constraint, Force, Velocity};
use crate::resources::*;
use bevy::prelude::*;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct PhysicsSet;

pub struct SimulationPlugin;

impl Plugin for SimulationPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(
            Update,
            (
                integrate_position_and_half_velocity,
                reset_forces,
                reset_energy,
                calculate_bond_and_constraint_forces,
                calculate_angle_forces,
                calculate_dihedral_forces,
                calculate_non_bonded_forces,
                sum_total_forces,
                finish_velocity_update,
                apply_thermostat,
                calculate_kinetic_energy,
                increment_step_count,
                end_simulation_step,
            )
                .chain()
                .in_set(PhysicsSet)
                .run_if(in_state(AppState::Running)),
        );
        app.add_systems(
            Update,
            (
                reset_forces,
                reset_energy,
                calculate_bond_and_constraint_forces,
                calculate_angle_forces,
                calculate_dihedral_forces,
                calculate_non_bonded_forces,
                sum_total_forces,
                log_initial_state_and_pause,
            )
                .chain()
                .run_if(in_state(AppState::PreSimulation)),
        );
    }
}

fn increment_step_count(mut step_count: ResMut<StepCount>) {
    step_count.0 += 1;
}

fn calculate_dihedral_forces(
    connectivity: Res<SystemConnectivity>,
    mut atom_query: Query<(&mut Transform, &mut Force, &Atom)>,
    force_field: Res<ForceField>,
    mut energy: ResMut<SystemEnergy>,
) {
    for dihedral in &connectivity.dihedrals {
        let Ok(
            [
                (t_a, mut f_a, type_a),
                (t_b, mut f_b, type_b),
                (t_c, mut f_c, type_c),
                (t_d, mut f_d, type_d),
            ],
        ) = atom_query.get_many_mut([dihedral.a, dihedral.b, dihedral.c, dihedral.d])
        else {
            continue;
        };

        let key = (
            type_a.type_name.clone(),
            type_b.type_name.clone(),
            type_c.type_name.clone(),
            type_d.type_name.clone(),
        );

        if let Some(&(k, n, phi0)) = force_field.dihedral_params.get(&key) {
            let r_ab = t_b.translation - t_a.translation;
            let r_cb = t_b.translation - t_c.translation;
            let r_dc = t_c.translation - t_d.translation;

            let n1 = r_ab.cross(r_cb).normalize_or_zero();
            let n2 = r_dc.cross(r_cb).normalize_or_zero();

            let phi = n1.angle_between(n2);
            let sign = if r_ab.dot(n2) > 0.0 { 1.0 } else { -1.0 };
            let phi = phi * sign;

            // Potential energy: E = k * (1 + cos(n*phi - phi0))
            energy.potential += k * (1.0 + (n as f32 * phi - phi0).cos());

            // Force calculation
            let force_mag = k * n as f32 * (n as f32 * phi - phi0).sin();

            let r_cb_len = r_cb.length();
            let f1 = (force_mag * r_cb_len / r_ab.length_squared()) * n1;
            let f4 = (-force_mag * r_cb_len / r_dc.length_squared()) * n2;

            let term_b = r_ab.dot(r_cb) / r_cb.length_squared();
            let term_c = r_dc.dot(r_cb) / r_cb.length_squared();

            let f2 = (1.0 - term_b) * f1 - term_c * f4;
            let f3 = (1.0 - term_c) * f4 - term_b * f1;

            f_a.angle += f1;
            f_b.angle += f2;
            f_c.angle += f3;
            f_d.angle += f4;
        }
    }
}

fn reset_forces(mut query: Query<&mut Force>) {
    for mut force in &mut query {
        *force = Force::default();
    }
}

fn apply_thermostat(
    mut query: Query<(&mut Velocity, &Atom)>,
    params: Res<SimulationParameters>,
    thermostat: Res<Thermostat>,
    force_field: Res<ForceField>,
    mut current_temp_res: ResMut<CurrentTemperature>,
    mut thermostat_scale_res: ResMut<ThermostatScale>,
) {
    let mut total_ke = 0.0;
    let mut num_atoms = query.iter().len();

    for (velocity, atom) in &query {
        // SAFE LOOKUP
        if let Some(atom_params) = force_field.atom_types.get(&atom.type_name) {
            total_ke += 0.5 * atom_params.mass * velocity.length_squared();
            num_atoms += 1;
        } else {
            warn!(
                "Skipping atom with unknown type '{}' in thermostat calc",
                atom.type_name
            );
        }
    }
    if num_atoms == 0 {
        return;
    }

    // Degrees of Freedom: 3 for each atom, minus 3 for the center of mass motion.
    let dof = (3 * num_atoms - 3) as f32;
    // The ideal gas constant R in units compatible with your simulation: kJ/(mol·K)
    const R: f32 = 0.0083144621;
    if total_ke < 1e-9 {
        return;
    }
    let current_temp = (2.0 * total_ke) / (dof * R);
    current_temp_res.0 = current_temp;
    if thermostat.target_temperature <= 0.0 {
        return;
    }
    if current_temp < 1e-9 {
        return;
    }
    let term_inside_sqrt =
        1.0 + (params.dt / thermostat.tau) * ((thermostat.target_temperature / current_temp) - 1.0);
    let scale = term_inside_sqrt.max(0.0).sqrt();
    thermostat_scale_res.0 = scale;
    for (mut velocity, _) in &mut query {
        velocity.0 *= scale;
    }
}

fn end_simulation_step(mut step_count: ResMut<StepCount>) {
    step_count.0 += 1;
    // info!("Advanced to step {}", step_count.0);
}

fn integrate_position_and_half_velocity(
    mut query: Query<(&mut Transform, &mut Velocity, &Acceleration)>,
    params: Res<SimulationParameters>,
) {
    for (mut transform, mut velocity, acceleration) in &mut query {
        velocity.0 += 0.5 * acceleration.0 * params.dt;
        transform.translation += velocity.0 * params.dt;
    }
}

fn finish_velocity_update(
    mut query: Query<(&mut Velocity, &mut Acceleration, &Force, &Atom)>,
    params: Res<SimulationParameters>,
    force_field: Res<ForceField>,
) {
    for (mut velocity, mut acceleration, force, atom) in &mut query {
        if let Some(atom_params) = force_field.atom_types.get(&atom.type_name) {
            let new_acceleration = force.total / atom_params.mass;
            velocity.0 += 0.5 * new_acceleration * params.dt;
            acceleration.0 = new_acceleration;
        } else {
            warn!(
                "Skipping atom with unknown type '{}' in velocity update",
                atom.type_name
            );
        }
    }
}

fn calculate_non_bonded_forces(
    mut query: Query<(Entity, &Atom, &GlobalTransform, &mut Force)>,
    force_field: Res<ForceField>,
    excluded: Res<ExcludedPairs>,
    mut energy: ResMut<SystemEnergy>,
) {
    let mut combinations = query.iter_combinations_mut();
    while let Some(
        [
            (e1, atom1_comp, transform1, mut force1),
            (e2, atom2_comp, transform2, mut force2),
        ],
    ) = combinations.fetch_next()
    {
        let key = if e1 < e2 { (e1, e2) } else { (e2, e1) };
        // Skip 1-2 pairs entirely
        if excluded.one_two.contains(&key) {
            continue;
        }
        let vec = transform2.translation() - transform1.translation();

        let r_sq = vec.length_squared();
        if r_sq < 1e-6 {
            continue;
        }
        let r = r_sq.sqrt();

        // Use atom-type specific parameters with Lorentz-Berthelot mixing rules
        let (Some(p1), Some(p2)) = (
            force_field.atom_types.get(&atom1_comp.type_name),
            force_field.atom_types.get(&atom2_comp.type_name),
        ) else {
            // If either atom type is unknown, we can't calculate non-bonded forces.
            continue;
        };
        let epsilon = (p1.epsilon * p2.epsilon).sqrt();
        let sigma = (p1.sigma + p2.sigma) * 0.5;

        // Only calculate LJ if epsilon is non-zero
        let lj_force_mag = if epsilon > 0.0 {
            let sigma_over_r = sigma / r;
            let sigma_over_r_6 = sigma_over_r.powi(6);
            let sigma_over_r_12 = sigma_over_r_6.powi(2);
            // LJ Potential: E = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]

            energy.potential += 4.0 * epsilon * (sigma_over_r_12 - sigma_over_r_6);
            (24.0 * epsilon / r) * (2.0 * sigma_over_r_12 - sigma_over_r_6)
        } else {
            0.0 // Skip LJ calculation if either atom has epsilon=0
        };

        // -- Electrostatics (Coulomb's Law) --
        let q1 = p1.charge;
        let q2 = p2.charge;

        let coulomb_scale = if excluded.one_three.contains(&key) {
            0.5 // Scaling factor for 1-3 pairs
        } else {
            1.0 // Full strength for other pairs
        };
        energy.potential += coulomb_scale * force_field.coulomb_k * (q1 * q2 / r);

        let coulomb_force_mag = force_field.coulomb_k * q1 * q2 / r_sq * coulomb_scale;
        let total_force_mag = lj_force_mag + coulomb_force_mag;
        let force_vec = vec.normalize_or_zero() * total_force_mag;

        force1.non_bonded -= force_vec;
        force2.non_bonded += force_vec;
    }
}

fn reset_energy(mut energy: ResMut<SystemEnergy>) {
    *energy = SystemEnergy::default();
}

fn calculate_kinetic_energy(
    mut energy: ResMut<SystemEnergy>,
    query: Query<(&Atom, &Velocity)>,
    force_field: Res<ForceField>,
    step_count: ResMut<StepCount>,
    mut history: ResMut<EnergyHistory>,
) {
    let mut kinetic_energy = 0.0;
    for (atom_type, velocity) in &query {
        if let Some(atom_params) = force_field.atom_types.get(&atom_type.type_name) {
            kinetic_energy += 0.5 * atom_params.mass * velocity.length_squared();
        } else {
            warn!(
                "Skipping atom with unknown type '{}' in kinetic energy calc",
                atom_type.type_name
            );
        }
    }
    energy.kinetic = kinetic_energy;
    energy.total = energy.potential + energy.kinetic;

    let step = step_count.0 as f64;
    history.potential.push_back((step, energy.potential as f64));
    history.kinetic.push_back((step, energy.kinetic as f64));
    history.total.push_back((step, energy.total as f64));

    if history.potential.len() > history.capacity {
        history.potential.pop_front();
        history.kinetic.pop_front();
        history.total.pop_front();
    }
}

fn sum_total_forces(mut query: Query<&mut Force>, multiplier: Res<ForceMultiplier>) {
    for mut force in &mut query {
        force.total = (force.bond + force.angle + force.non_bonded) * multiplier.0;
    }
}

// System 3: Calculate angle forces
fn calculate_angle_forces(
    connectivity: Res<SystemConnectivity>,
    mut atom_query: Query<(&Transform, &mut Force, &Atom)>,
    force_field: Res<ForceField>,
    mut energy: ResMut<SystemEnergy>,
) {
    for angle in &connectivity.angles {
        let Ok(
            [
                (t1, mut f1, type1),
                (t_center, mut f_center, type_center),
                (t2, mut f2, type2),
            ],
        ) = atom_query.get_many_mut([angle.a, angle.center, angle.b])
        else {
            continue;
        };

        if let Some(&(angle_k, angle_theta0)) = force_field.angle_params.get(&(
            type1.type_name.clone(),
            type_center.type_name.clone(),
            type2.type_name.clone(),
        )) {
            let v1 = t1.translation - t_center.translation;
            let v2 = t2.translation - t_center.translation;
            let theta = v1.angle_between(v2);

            let force_magnitude = -angle_k * (theta - angle_theta0);

            let force1_dir = v1.cross(v1.cross(v2)).normalize_or_zero();
            let force2_dir = v2.cross(v2.cross(v1)).normalize_or_zero();

            let force_on_1 = force1_dir * force_magnitude;
            let force_on_2 = force2_dir * force_magnitude;
            f1.angle += force_on_1;
            f2.angle += force_on_2;
            f_center.angle -= force_on_1 + force_on_2;

            energy.potential += 0.5 * angle_k * (theta - angle_theta0).powi(2);
        }
    }
}

fn log_initial_state_and_pause(
    mut next_state: ResMut<NextState<AppState>>,
    atom_query: Query<(Entity, &Transform, &Velocity, &Force), With<Atom>>,
) {
    info!("--- INITIAL STATE LOG ---");
    info!("Atom Count: {}", atom_query.iter().len());
    for (entity, transform, velocity, force) in &atom_query {
        info!(
            "  Entity: {:?}, Pos: {:?}, Vel: {:?}, Force: {:?}",
            entity, transform.translation, velocity.0, force.total
        );

        // Check for non-finite numbers which cause explosions
        if !transform.translation.is_finite() {
            error!("!!! NON-FINITE POSITION on Entity {:?}", entity);
        }
        if !velocity.0.is_finite() {
            error!("!!! NON-FINITE VELOCITY on Entity {:?}", entity);
        }
        if !force.total.is_finite() {
            error!("!!! NON-FINITE FORCE on Entity {:?}", entity);
        }
    }
    info!("--- END INITIAL STATE LOG ---");
    info!("Pausing simulation for inspection.");

    // Force the app to pause so we can read the logs.
    next_state.set(AppState::Paused);
}

fn calculate_bond_and_constraint_forces(
    connectivity: Res<SystemConnectivity>,
    drag_state: Res<DragState>, // Just need to read the state
    constraint_q: Query<(Entity, &Constraint)>,
    mut atom_query: Query<(&Transform, &mut Force, &Atom)>,
    force_field: Res<ForceField>,
    mut energy: ResMut<SystemEnergy>,
) {
    // --- 1. The original bond force calculation loop (unchanged) ---
    for bond in &connectivity.bonds {
        if let Ok([(t1, mut f1, type1), (t2, mut f2, type2)]) =
            atom_query.get_many_mut([bond.a, bond.b])
        {
            if let Some(&(bond_k, bond_r0)) = force_field.bond_params.get(&(
                type1.type_name.clone(),
                type2.type_name.clone(),
                bond.order,
            )) {
                let vec = t2.translation - t1.translation;

                let r = vec.length();
                let force_magnitude = -bond_k * (r - bond_r0);
                let force_vec = vec.normalize_or_zero() * force_magnitude;
                f1.bond -= force_vec;
                f2.bond += force_vec;
                energy.potential += 0.5 * bond_k * (r - bond_r0).powi(2);
            }
        }
    }

    if let (Some(initial_centroid), Some(target_centroid)) =
        (drag_state.initial_centroid, drag_state.target_centroid)
    {
        let delta = target_centroid - initial_centroid;

        for (i, (entity, constraint)) in constraint_q.iter().enumerate() {
            if let Ok((atom_transform, mut force, _)) = atom_query.get_mut(entity) {
                let target_pos = drag_state.initial_positions[i] + delta;
                let force_vec = (target_pos - atom_transform.translation) * constraint.stiffness;
                force.bond += force_vec;
            }
        }
    }
}
