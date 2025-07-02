use crate::components::{Atom, BondVisualization};
use bevy::prelude::*;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct VisualizationSet;

pub struct VisualizationPlugin;

impl Plugin for VisualizationPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Update, update_bond_visuals.in_set(VisualizationSet));
    }
}

fn update_bond_visuals(
    // Query for all atom transforms (read-only)
    atom_query: Query<&Transform, With<Atom>>,
    // Query for the bond entities we need to update
    mut bond_query: Query<(&mut Transform, &BondVisualization), Without<Atom>>,
) {
    const BOND_SPACING: f32 = 0.035;
    for (mut bond_transform, bond_vis) in &mut bond_query {
        // Get the positions of the two atoms for this bond
        if let Ok([t1, t2]) = atom_query.get_many([bond_vis.atom1, bond_vis.atom2]) {
            let start_pos = t1.translation;
            let end_pos = t2.translation;

            let distance = start_pos.distance(end_pos);
            bond_transform.scale = Vec3::new(1.0, distance, 1.0);
            let direction = (end_pos - start_pos).normalize_or_zero();
            bond_transform.rotation = Quat::from_rotation_arc(Vec3::Y, direction);

            // Offset logic for multi-bond
            let midpoint = (start_pos + end_pos) / 2.0;
            if bond_vis.total_strands <= 1 {
                bond_transform.translation = midpoint;
            } else {
                // 1. Find a vector perpendicular to the bond's direction.
                // We use a common trick to avoid the case where the bond is aligned with our 'up' vector.
                let up_vector = if direction.dot(Vec3::Y).abs() < 0.999 {
                    Vec3::Z
                } else {
                    Vec3::X
                };
                let offset_dir = direction.cross(up_vector).normalize_or_zero();

                // 2. Calculate how far this specific strand should be from the center.
                let total_width = BOND_SPACING * (bond_vis.total_strands - 1) as f32;
                let strand_offset_dist =
                    BOND_SPACING * bond_vis.strand_index as f32 - total_width / 2.0;

                // 3. Apply the final offset.
                bond_transform.translation = midpoint + (offset_dir * strand_offset_dist);
            }
        }
    }
}
