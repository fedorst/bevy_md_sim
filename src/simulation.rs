use super::core::*;
use bevy::prelude::*;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct PhysicsSet;

pub struct SimulationPlugin;

impl Plugin for SimulationPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(
            Update,
            (
                (
                    integrate_position_and_half_velocity,
                    reset_forces,
                    calculate_bond_forces,
                    calculate_angle_forces,
                    calculate_non_bonded_forces,
                    finish_velocity_update,
                    end_simulation_step,
                )
                    .chain()
                    .in_set(PhysicsSet)
                    .run_if(
                        |state: Res<SimulationState>, step: Res<StepSimulation>| {
                            !state.paused || step.0
                        },
                    ),
                end_simulation_step.after(PhysicsSet),
            ),
        );
    }
}
fn reset_forces(mut query: Query<&mut Force>) {
    for mut force in &mut query {
        force.0 = Vec3::ZERO;
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
        // Half-step velocity update
        velocity.0 += 0.5 * acceleration.0 * params.dt;
        // Full-step position update
        transform.translation += velocity.0 * params.dt;
    }
}

fn finish_velocity_update(
    mut query: Query<(&mut Velocity, &mut Acceleration, &Force, &AtomType)>,
    params: Res<SimulationParameters>,
) {
    for (mut velocity, mut acceleration, force, atom_type) in &mut query {
        // Calculate and store the new acceleration
        let new_acceleration = force.0 / atom_type.mass();
        // Complete the velocity update
        velocity.0 += 0.5 * new_acceleration * params.dt;
        // Store the new acceleration for the next step's position update
        acceleration.0 = new_acceleration;
    }
}

fn calculate_non_bonded_forces(
    mut query: Query<(Entity, &AtomType, &GlobalTransform, &mut Force)>,
    force_field: Res<ForceField>,
    excluded: Res<ExcludedPairs>,
) {
    let mut combinations = query.iter_combinations_mut();
    while let Some(
        [
            (e1, atom1, transform1, mut force1),
            (e2, atom2, transform2, mut force2),
        ],
    ) = combinations.fetch_next()
    {
        // ignore non-bonded for 1-2 and 1-3
        let key = if e1 < e2 { (e1, e2) } else { (e2, e1) };
        if excluded.0.contains(&key) {
            // info!("Skipping non-bonded for {:?} and {:?}", e1, e2);
            continue;
        }
        let vec = transform2.translation() - transform1.translation();
        let r_sq = vec.length_squared();
        if r_sq < 0.0001 {
            continue;
        } // Avoid division by zero
        let r = r_sq.sqrt();

        /*
        info!(
            "Non-bonded {:?}({:?}) - {:?}({:?}): r={:.6}",
            e1, atom1, e2, atom2, r
        );
        */

        // Use atom-type specific parameters with Lorentz-Berthelot mixing rules
        let epsilon = (atom1.epsilon() * atom2.epsilon()).sqrt();
        let sigma = (atom1.sigma() + atom2.sigma()) * 0.5;

        // Only calculate LJ if epsilon is non-zero
        let lj_force_mag = if epsilon > 0.0 {
            let sigma_over_r = sigma / r;
            let sigma_over_r_6 = sigma_over_r.powi(6);
            let sigma_over_r_12 = sigma_over_r_6.powi(2);
            (24.0 * epsilon / r) * (2.0 * sigma_over_r_12 - sigma_over_r_6)
        } else {
            0.0 // Skip LJ calculation if either atom has epsilon=0
        };

        // -- Electrostatics (Coulomb's Law) --
        let q1 = atom1.charge();
        let q2 = atom2.charge();
        let coulomb_force_mag = force_field.coulomb_k * q1 * q2 / r_sq;

        let total_force_mag = lj_force_mag + coulomb_force_mag;
        let force_vec = vec.normalize_or_zero() * total_force_mag;

        force1.0 -= force_vec;
        force2.0 += force_vec;

        /*
        info!(
            "  LJ: F={:.6}, Coulomb: F={:.6}",
            lj_force_mag, coulomb_force_mag
        );
         */
    }
}

// System 3: Calculate angle forces
fn calculate_angle_forces(
    connectivity: Res<SystemConnectivity>,
    mut atom_query: Query<(&Transform, &mut Force, &AtomType)>,
    force_field: Res<ForceField>,
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
            f1.0 += force_on_1;
            f2.0 += force_on_2;
            f_center.0 -= force_on_1 + force_on_2;
        }

        /*
        info!(
            "Angle {:?}-{:?}-{:?}: θ={:.2}°, Δθ={:.2}°, F={:.6}",
            angle.a,
            angle.center,
            angle.b,
            theta.to_degrees(),
            (theta - force_field.angle_theta0).to_degrees(),
            force_magnitude
        );  */
    }
}

fn calculate_bond_forces(
    connectivity: Res<SystemConnectivity>,
    mut atom_query: Query<(&Transform, &mut Force, &AtomType)>,
    force_field: Res<ForceField>,
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
            f1.0 -= force_vec;
            f2.0 += force_vec;
            /*
            info!(
                "Bond {:?}-{:?}: r={:.6}, Δr={:.6}, F={:.6}",
                bond.a,
                bond.b,
                r,
                r - force_field.bond_r0,
                force_magnitude
            );
             */
        }
    }
}
