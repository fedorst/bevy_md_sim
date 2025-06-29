use super::core::*;
use bevy::prelude::*;

pub struct VisualizationPlugin;

impl Plugin for VisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Update, update_bond_visuals);
    }
}

fn update_bond_visuals(
    // Query for all atom transforms (read-only)
    atom_query: Query<&Transform, With<AtomType>>,
    // Query for the bond entities we need to update
    mut bond_query: Query<(&mut Transform, &BondVisualization), Without<AtomType>>,
) {
    for (mut bond_transform, bond_vis) in &mut bond_query {
        // Get the positions of the two atoms for this bond
        if let (Ok(t1), Ok(t2)) = (
            atom_query.get(bond_vis.atom1),
            atom_query.get(bond_vis.atom2),
        ) {
            let start_pos = t1.translation;
            let end_pos = t2.translation;

            // 1. Position: The center of the cylinder is the midpoint
            let midpoint = (start_pos + end_pos) / 2.0;
            bond_transform.translation = midpoint;

            // 2. Scale: The length of the cylinder is the distance between atoms
            let distance = start_pos.distance(end_pos);
            // The cylinder mesh has a height of 1.0, so we scale it by the distance.
            bond_transform.scale = Vec3::new(1.0, distance, 1.0);

            // 3. Rotation: Make the cylinder's local Y-axis point from start to end
            let direction = (end_pos - start_pos).normalize_or_zero();
            // `Quat::from_rotation_arc` finds the rotation between two vectors.
            // We want to rotate the cylinder's default up-axis (Vec3::Y) to our new direction.
            bond_transform.rotation = Quat::from_rotation_arc(Vec3::Y, direction);
        }
    }
}
