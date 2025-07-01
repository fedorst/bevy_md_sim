use crate::components::{Acceleration, AtomType, Force, Velocity};
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
                calculate_bond_forces,
                calculate_angle_forces,
                calculate_non_bonded_forces,
                sum_total_forces,
                finish_velocity_update,
                apply_thermostat,
                calculate_kinetic_energy,
                end_simulation_step,
            )
                .chain()
                .in_set(PhysicsSet)
                .run_if(|state: Res<SimulationState>, step: Res<StepSimulation>| {
                    !state.paused || step.0
                }),
        );
    }
}
fn reset_forces(mut query: Query<&mut Force>) {
    for mut force in &mut query {
        *force = Force::default();
    }
}

fn apply_thermostat(
    mut query: Query<(&mut Velocity, &AtomType)>,
    params: Res<SimulationParameters>,
    thermostat: Res<Thermostat>,
    mut current_temp_res: ResMut<CurrentTemperature>,
    mut thermostat_scale_res: ResMut<ThermostatScale>,
) {
    let mut total_ke = 0.0;
    let num_atoms = query.iter().len();
    if num_atoms == 0 {
        return;
    }
    for (velocity, atom_type) in &query {
        total_ke += 0.5 * atom_type.mass() * velocity.length_squared();
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

fn end_simulation_step(mut step: ResMut<StepSimulation>, mut step_count: ResMut<StepCount>) {
    step_count.0 += 1;
    step.0 = false;
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
    mut query: Query<(&mut Velocity, &mut Acceleration, &Force, &AtomType)>,
    params: Res<SimulationParameters>,
) {
    for (mut velocity, mut acceleration, force, atom_type) in &mut query {
        let new_acceleration = force.total / atom_type.mass();
        velocity.0 += 0.5 * new_acceleration * params.dt;
        acceleration.0 = new_acceleration;
    }
}

fn calculate_non_bonded_forces(
    mut query: Query<(Entity, &AtomType, &GlobalTransform, &mut Force)>,
    force_field: Res<ForceField>,
    excluded: Res<ExcludedPairs>,
    mut energy: ResMut<SystemEnergy>,
) {
    let mut combinations = query.iter_combinations_mut();
    while let Some(
        [
            (e1, atom1, transform1, mut force1),
            (e2, atom2, transform2, mut force2),
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
        let epsilon = (atom1.epsilon() * atom2.epsilon()).sqrt();
        let sigma = (atom1.sigma() + atom2.sigma()) * 0.5;

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
        let q1 = atom1.charge();
        let q2 = atom2.charge();

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
    query: Query<(&AtomType, &Velocity)>,
) {
    let mut kinetic_energy = 0.0;
    for (atom_type, velocity) in &query {
        kinetic_energy += 0.5 * atom_type.mass() * velocity.length_squared();
    }
    energy.kinetic = kinetic_energy;
    energy.total = energy.potential + energy.kinetic;
}

fn sum_total_forces(mut query: Query<&mut Force>) {
    for mut force in &mut query {
        force.total = force.bond + force.angle + force.non_bonded;
    }
}

// System 3: Calculate angle forces
fn calculate_angle_forces(
    connectivity: Res<SystemConnectivity>,
    mut atom_query: Query<(&Transform, &mut Force, &AtomType)>,
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

        if let Some(&(angle_k, angle_theta0)) =
            force_field
                .angle_params
                .get(&(*type1, *type_center, *type2))
        {
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

fn calculate_bond_forces(
    connectivity: Res<SystemConnectivity>,
    mut atom_query: Query<(&Transform, &mut Force, &AtomType)>,
    force_field: Res<ForceField>,
    mut energy: ResMut<SystemEnergy>,
) {
    for bond in &connectivity.bonds {
        let Ok([(t1, mut f1, type1), (t2, mut f2, type2)]) =
            atom_query.get_many_mut([bond.a, bond.b])
        else {
            continue;
        };

        if let Some(&(bond_k, bond_r0)) = force_field.bond_params.get(&(*type1, *type2)) {
            let vec = t2.translation - t1.translation;
            let r = vec.length();
            let force_magnitude = -bond_k * (r - bond_r0);
            let force_vec = vec.normalize_or_zero() * force_magnitude;
            f1.bond -= force_vec;
            f2.bond += force_vec;
            // for accum
            energy.potential += 0.5 * bond_k * (r - bond_r0).powi(2);
        }
    }
}
