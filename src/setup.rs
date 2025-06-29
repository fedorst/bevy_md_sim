use super::core::*;
use bevy::prelude::*;
use bevy::ui::{PositionType, Val};
use std::collections::HashSet;
#[derive(Component)]
pub struct PauseText;

use bevy_panorbit_camera::PanOrbitCamera;

pub struct SetupPlugin;

impl Plugin for SetupPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, (setup, build_exclusions).chain());
    }
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut connectivity: ResMut<SystemConnectivity>,
) {
    commands.spawn((
        PanOrbitCamera::default(),
        Transform::from_xyz(-1.0, 2.0, 2.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        DirectionalLight::default(),
        Transform::from_xyz(10.0, 20.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    commands.spawn((
        Text::new("PAUSED"), // Spawn a Text component directly
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(10.0),
            right: Val::Px(10.0),
            ..default()
        },
        TextFont {
            font_size: 30.0,
            ..default()
        }, // Add font component
        TextColor(Color::WHITE), // Add color component
        PauseText,               // Add our marker
    ));

    let c_mesh = meshes.add(Sphere::new(0.07));
    let o_mesh = meshes.add(Sphere::new(0.06));
    let h_mesh = meshes.add(Sphere::new(0.04));

    let atom_data = [
        // AtomType, Position (nm)
        (AtomType::Carbon, Vec3::new(-0.0658, 0.038, 0.0)), // C1
        (AtomType::Carbon, Vec3::new(0.0658, -0.038, 0.0)), // C2
        (AtomType::Oxygen, Vec3::new(0.133, -0.05, 0.1224)), // O
        (AtomType::Hydrogen, Vec3::new(-0.124, 0.09, 0.089)), // H on C1
        (AtomType::Hydrogen, Vec3::new(-0.124, 0.09, -0.089)), // H on C1
        (AtomType::Hydrogen, Vec3::new(-0.074, -0.06, 0.0)), // H on C1
        (AtomType::Hydrogen, Vec3::new(0.074, 0.06, 0.089)), // H on C2
        (AtomType::Hydrogen, Vec3::new(0.074, 0.06, -0.089)), // H on C2
        (AtomType::Hydrogen, Vec3::new(0.18, -0.1, 0.04)),  // H on O
    ];

    let mut atoms = Vec::new();
    for (atom_type, pos) in atom_data {
        let (mesh, mat) = match atom_type {
            AtomType::Carbon => (c_mesh.clone(), materials.add(Color::srgb(0.2, 0.2, 0.2))),
            AtomType::Oxygen => (o_mesh.clone(), materials.add(Color::srgb(1.0, 0.1, 0.1))),
            AtomType::Hydrogen => (h_mesh.clone(), materials.add(Color::srgb(0.9, 0.9, 0.9))),
        };
        let entity = commands
            .spawn((
                atom_type,
                Force::default(),
                Velocity(Vec3::ZERO),
                Acceleration(Vec3::ZERO),
                Mesh3d(mesh),
                MeshMaterial3d(mat),
                Transform::from_translation(pos),
            ))
            .id();
        atoms.push(entity);
    }

    let bonds_to_add = [
        (0, 1),
        (0, 3),
        (0, 4),
        (0, 5), // C1 bonds
        (1, 2),
        (1, 6),
        (1, 7), // C2 bonds
        (2, 8), // O-H bond
    ];
    for (i, j) in bonds_to_add {
        connectivity.bonds.push(Bond {
            a: atoms[i],
            b: atoms[j],
        });
    }

    let angles_to_add = [
        (1, 0, 3),
        (1, 0, 4),
        (1, 0, 5),
        (3, 0, 4),
        (3, 0, 5),
        (4, 0, 5), // Angles around C1
        (0, 1, 2),
        (0, 1, 6),
        (0, 1, 7),
        (2, 1, 6),
        (2, 1, 7),
        (6, 1, 7), // Angles around C2
        (1, 2, 8), // C-O-H angle
    ];
    for (i, j, k) in angles_to_add {
        connectivity.angles.push(Angle {
            a: atoms[i],
            center: atoms[j],
            b: atoms[k],
        });
    }
    let bond_mesh = meshes.add(Cylinder::new(0.02, 1.0)); // Radius 0.02nm, length will be set by scale
    let bond_mat = materials.add(Color::srgb(0.8, 0.8, 0.2)); // Yellow

    for bond in &connectivity.bonds {
        commands.spawn((
            Mesh3d(bond_mesh.clone()),
            MeshMaterial3d(bond_mat.clone()),
            // default Transform. The update system will set its value each frame.
            Transform::default(),
            BondVisualization {
                atom1: bond.a,
                atom2: bond.b,
            },
        ));
    }
}

fn build_exclusions(mut commands: Commands, connectivity: Res<SystemConnectivity>) {
    let mut excluded = HashSet::new();
    // filter out direct bonds from non-bonded interactions
    for bond in &connectivity.bonds {
        let key = if bond.a < bond.b {
            (bond.a, bond.b)
        } else {
            (bond.b, bond.a)
        };
        excluded.insert(key);
    }

    // THIS IS THE FIX: Correctly add all 1-3 pairs (from angles)
    // The 1-3 pair in an angle A-C-B are the two outer atoms, A and B.
    for angle in &connectivity.angles {
        let key = if angle.a < angle.b {
            (angle.a, angle.b)
        } else {
            (angle.b, angle.a)
        };
        excluded.insert(key);
    }
    info!("Built exclusion list with {} pairs.", excluded.len());
    commands.insert_resource(ExcludedPairs(excluded));
}
