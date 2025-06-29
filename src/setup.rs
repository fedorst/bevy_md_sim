use super::config::MoleculeConfig;
use super::core::*;
use bevy::prelude::*;
use bevy::ui::{PositionType, Val};
use bevy_panorbit_camera::PanOrbitCamera;
use std::collections::{HashMap, HashSet};

#[derive(Component)]
pub struct PauseText;

#[derive(Component)]
pub struct TempDisplayText;

#[derive(Component)]
pub struct TimeDisplayText;

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
    molecule_selection: Res<MoleculeSelection>,
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

    info!("Loading molecule from config file...");
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let filename = format!("{}.json", molecule_selection.0);
    let config_path = std::path::Path::new(&manifest_dir)
        .join("assets/molecules")
        .join(filename);

    let file_contents = std::fs::read_to_string(&config_path).expect(&format!(
        "Failed to read molecule file at {:?}",
        config_path
    ));
    let molecule_config: MoleculeConfig =
        serde_json::from_str(&file_contents).expect("Failed to parse molecule JSON");

    info!("Successfully loaded molecule: {}", molecule_config.name);

    let bond_mesh = meshes.add(Cylinder::new(0.02, 1.0));
    let bond_mat = materials.add(Color::srgb(0.8, 0.8, 0.2));

    let atom_mesh_map = [
        (AtomType::Carbon, meshes.add(Sphere::new(0.07))),
        (AtomType::Oxygen, meshes.add(Sphere::new(0.06))),
        (AtomType::Hydrogen, meshes.add(Sphere::new(0.04))),
        (AtomType::Nitrogen, meshes.add(Sphere::new(0.065))),
    ]
    .into_iter()
    .collect::<HashMap<_, _>>();

    let mut id_to_entity_map: HashMap<String, Entity> = HashMap::new();

    for atom_spec in &molecule_config.atoms {
        let core_atom_type: AtomType = atom_spec.atom_type.into();

        let entity = commands
            .spawn((
                core_atom_type,
                Force::default(),
                Velocity(Vec3::ZERO),
                Acceleration(Vec3::ZERO),
                Mesh3d(atom_mesh_map[&core_atom_type].clone()),
                MeshMaterial3d(materials.add(match core_atom_type {
                    AtomType::Carbon => Color::srgb(0.2, 0.2, 0.2),
                    AtomType::Oxygen => Color::srgb(1.0, 0.1, 0.1),
                    AtomType::Hydrogen => Color::srgb(0.9, 0.9, 0.9),
                    AtomType::Nitrogen => Color::srgb(0.1, 0.1, 1.0),
                })),
                Transform::from_translation(Vec3::from(atom_spec.pos)),
            ))
            .id();

        // Populate our ID-to-Entity map.
        id_to_entity_map.insert(atom_spec.id.clone(), entity);
    }

    // --- 2. Create bonds from the config file ---
    for bond_spec in &molecule_config.bonds {
        let entity1 = id_to_entity_map
            .get(&bond_spec[0])
            .expect("Atom ID not found in map");
        let entity2 = id_to_entity_map
            .get(&bond_spec[1])
            .expect("Atom ID not found in map");

        connectivity.bonds.push(Bond {
            a: *entity1,
            b: *entity2,
        });
    }

    // --- 3. ALGORITHMICALLY Generate Angles and Bond Visuals ---
    info!("Generating angles and bond visuals...");

    let mut adjacency: HashMap<Entity, Vec<Entity>> = HashMap::new();
    for bond in &connectivity.bonds {
        adjacency.entry(bond.a).or_default().push(bond.b);
        adjacency.entry(bond.b).or_default().push(bond.a);

        // Spawn the visual cylinder for the bond
        commands.spawn((
            Mesh3d(bond_mesh.clone()),
            MeshMaterial3d(bond_mat.clone()),
            Transform::default(),
            BondVisualization {
                atom1: bond.a,
                atom2: bond.b,
            },
        ));
    }

    for (center_atom, neighbors) in adjacency {
        if neighbors.len() < 2 {
            continue;
        }
        // Generate all unique pairs of neighbors for the angles
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                connectivity.angles.push(Angle {
                    a: neighbors[i],
                    center: center_atom,
                    b: neighbors[j],
                });
            }
        }
    }
    info!(
        "Setup complete. Atoms: {}, Bonds: {}, Angles: {}",
        molecule_config.atoms.len(),
        connectivity.bonds.len(),
        connectivity.angles.len()
    );
}

fn build_exclusions(mut commands: Commands, connectivity: Res<SystemConnectivity>) {
    let mut one_two = HashSet::new();
    let mut one_three = HashSet::new();

    for bond in &connectivity.bonds {
        let key = if bond.a < bond.b {
            (bond.a, bond.b)
        } else {
            (bond.b, bond.a)
        };
        one_two.insert(key);
    }

    for angle in &connectivity.angles {
        let key = if angle.a < angle.b {
            (angle.a, angle.b)
        } else {
            (angle.b, angle.a)
        };
        one_three.insert(key);
    }
    info!(
        "Built exclusion lists: {} 1-2 pairs, {} 1-3 pairs",
        one_two.len(),
        one_three.len()
    );
    commands.insert_resource(ExcludedPairs { one_two, one_three });
}
