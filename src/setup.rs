use super::components::*;
use super::config::MoleculeConfig;
use super::resources::{BondOrder, *};
use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCamera;
use std::collections::{HashMap, HashSet};

pub struct SetupPlugin;

impl Plugin for SetupPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(
            Startup,
            (
                setup_scene,
                load_molecule_from_config,
                build_derived_connectivity,
                spawn_bond_visuals,
            )
                .chain(),
        );
    }
}

fn setup_scene(mut commands: Commands) {
    info!("Setting up scene: Camera, Light, UI");
    commands.spawn((
        PanOrbitCamera::default(),
        Transform::from_xyz(-1.0, 2.0, 2.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        DirectionalLight::default(),
        Transform::from_xyz(10.0, 20.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
}

fn load_molecule_from_config(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut connectivity: ResMut<SystemConnectivity>,
    force_field: Res<ForceField>,
    molecule_selection: Res<MoleculeSelection>,
    mut atom_id_map: ResMut<AtomIdMap>,
) {
    info!("Loading molecule from config file...");
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let filename = format!("{}.json", molecule_selection.0);
    let config_path = std::path::Path::new(&manifest_dir)
        .join("assets/molecules")
        .join(filename);

    let file_contents = std::fs::read_to_string(&config_path)
        .unwrap_or_else(|_| panic!("Failed to read molecule file at {:?}", config_path));
    let molecule_config: MoleculeConfig =
        serde_json::from_str(&file_contents).expect("Failed to parse molecule JSON");

    info!("Spawning atoms for: {}", molecule_config.name);

    let mut element_mesh_map: HashMap<String, Handle<Mesh>> = HashMap::new();

    let mut id_to_entity_map = HashMap::new();
    for atom_spec in &molecule_config.atoms {
        // let core_atom_type: AtomType = atom_spec.atom_type.into();
        let atom_params = force_field.atom_types.get(&atom_spec.type_name).unwrap();
        let element_str = &atom_params.element;

        let mesh_handle = element_mesh_map
            .entry(element_str.clone())
            .or_insert_with(|| {
                let radius = match element_str.as_str() {
                    "C" => 0.06,
                    "O" => 0.05,
                    "H" => 0.03,
                    "N" => 0.055,
                    _ => 0.045,
                };
                meshes.add(Sphere::new(radius))
            });

        let entity = commands
            .spawn((
                Atom {
                    type_name: atom_spec.type_name.clone(),
                },
                Force::default(),
                Velocity(Vec3::ZERO),
                Acceleration(Vec3::ZERO),
                Mesh3d(mesh_handle.clone()),
                MeshMaterial3d(materials.add(match element_str.as_str() {
                    "C" => Color::srgb(0.2, 0.2, 0.2),
                    "O" => Color::srgb(1.0, 0.1, 0.1),
                    "H" => Color::srgb(0.9, 0.9, 0.9),
                    "N" => Color::srgb(0.1, 0.1, 1.0),
                    _ => Color::srgb(1.0, 0.2, 0.8), // Default pink
                })),
                Transform::from_translation(Vec3::from(atom_spec.pos)),
            ))
            .id();

        id_to_entity_map.insert(atom_spec.id.clone(), entity);
        atom_id_map
            .entity_to_id
            .insert(entity, atom_spec.id.clone());
    }

    info!("Building bond list...");
    for bond_spec in &molecule_config.bonds {
        let entity1 = id_to_entity_map.get(&bond_spec[0]).unwrap();
        let entity2 = id_to_entity_map.get(&bond_spec[1]).unwrap();
        connectivity.bonds.push(Bond {
            a: *entity1,
            b: *entity2,
            order: BondOrder::Single,
        });
    }

    info!("Building double bond list (if any)...");
    if let Some(double_bonds) = molecule_config.double_bonds {
        for bond_spec in &double_bonds {
            let entity1 = id_to_entity_map.get(&bond_spec[0]).unwrap();
            let entity2 = id_to_entity_map.get(&bond_spec[1]).unwrap();
            connectivity.bonds.push(Bond {
                a: *entity1,
                b: *entity2,
                order: BondOrder::Double,
            });
        }
    }

    let bond_mesh_handle = meshes.add(Cylinder::new(0.015, 1.0));
    let bond_material_handle = materials.add(Color::srgb(0.8, 0.8, 0.2));

    commands.insert_resource(SharedAssetHandles {
        bond_mesh: bond_mesh_handle,
        bond_material: bond_material_handle,
    });
}

fn build_derived_connectivity(
    mut commands: Commands,
    mut connectivity: ResMut<SystemConnectivity>,
) {
    info!("Generating angles and exclusion lists...");
    let mut adjacency: HashMap<Entity, Vec<Entity>> = HashMap::new();
    for bond in &connectivity.bonds {
        adjacency.entry(bond.a).or_default().push(bond.b);
        adjacency.entry(bond.b).or_default().push(bond.a);
    }

    for (center_atom, neighbors) in adjacency {
        if neighbors.len() < 2 {
            continue;
        }
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

    let mut new_dihedrals = Vec::new();
    for bond1 in &connectivity.bonds {
        for bond2 in &connectivity.bonds {
            // Find dihedrals A-B-C-D where B-C is the connection between bond1 and bond2
            if bond1.b == bond2.a && bond1.a != bond2.b {
                new_dihedrals.push(Dihedral {
                    a: bond1.a,
                    b: bond1.b,
                    c: bond2.a,
                    d: bond2.b,
                });
            }
        }
    }
    // After the immutable borrows are finished, append the new dihedrals.
    connectivity.dihedrals.append(&mut new_dihedrals);

    // Now build the exclusion lists (this is the logic from your old `build_exclusions`)
    let mut one_two = HashSet::new();
    for bond in &connectivity.bonds {
        let key = if bond.a < bond.b {
            (bond.a, bond.b)
        } else {
            (bond.b, bond.a)
        };
        one_two.insert(key);
    }

    let mut one_three = HashSet::new();
    for angle in &connectivity.angles {
        let key = if angle.a < angle.b {
            (angle.a, angle.b)
        } else {
            (angle.b, angle.a)
        };
        one_three.insert(key);
    }

    info!(
        "Connectivity generation complete. Angles: {}, Dihedrals: {}, 1-2 Exclusions: {}, 1-3 Exclusions: {}",
        connectivity.angles.len(),
        connectivity.dihedrals.len(),
        one_two.len(),
        one_three.len()
    );
    commands.insert_resource(ExcludedPairs { one_two, one_three });
}

fn spawn_bond_visuals(
    mut commands: Commands,
    shared_handles: Res<SharedAssetHandles>,
    connectivity: Res<SystemConnectivity>,
) {
    info!("Spawning initial bond visuals...");
    for bond in &connectivity.bonds {
        let num_strands = match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
        };

        for i in 0..num_strands {
            commands.spawn((
                Mesh3d(shared_handles.bond_mesh.clone()),
                MeshMaterial3d(shared_handles.bond_material.clone()),
                Transform::default(),
                BondVisualization {
                    atom1: bond.a,
                    atom2: bond.b,
                    strand_index: i,
                    total_strands: num_strands,
                },
            ));
        }
    }
}
