// src/spawning.rs

use crate::components::*;
use crate::config::MoleculeConfig;
use crate::interaction::RebuildConnectivityEvent;
use crate::resources::{
    Angle, AtomIdMap, Bond, BondOrder, Dihedral, ExcludedPairs, ForceField, SharedAssetHandles,
    SystemConnectivity,
};
use bevy::prelude::*;
use bevy::tasks::{AsyncComputeTaskPool, Task};
use futures_lite::future;
use std::collections::{HashMap, HashSet};

/// An event sent to request a new molecule from a SMILES string.
#[derive(Event, Debug)]
pub struct SpawnMoleculeFromSMILESEvent(pub String);

#[derive(Event, Debug)]
pub struct ValidateSMILESEvent(pub String);
/// Event sent by the async task with the validation result.
#[derive(Event, Debug)]
pub struct SMILESValidationResult(pub Result<(), String>);

/// An event sent when the python script succeeds, containing the molecule JSON.
#[derive(Event, Debug)]
struct SpawnMoleculeFromJsonEvent(String);

/// A component to mark the container entity for all atoms and bonds of a molecule.
#[derive(Component)]
pub struct MoleculeContainer;

#[derive(Component)]
struct ValidationTask(Option<Task<SMILESValidationResult>>);

pub struct SpawningPlugin;

impl Plugin for SpawningPlugin {
    fn build(&self, app: &mut App) {
        app.add_event::<SpawnMoleculeFromSMILESEvent>()
            .add_event::<ValidateSMILESEvent>()
            .add_event::<SMILESValidationResult>()
            .add_event::<SpawnMoleculeFromJsonEvent>()
            .add_systems(
                Update,
                (
                    trigger_smiles_validation,
                    handle_validation_result,
                    trigger_molecule_generation_on_enter,
                    (
                        despawn_previous_molecule,
                        ApplyDeferred,
                        spawn_new_molecule,
                        build_derived_connectivity,
                    )
                        .chain()
                        .run_if(on_event::<SpawnMoleculeFromJsonEvent>),
                ),
            );
    }
}

fn trigger_smiles_validation(
    mut commands: Commands,
    mut validation_events: EventReader<ValidateSMILESEvent>,
    mut existing_task_query: Query<(Entity, &mut ValidationTask)>,
) {
    // Read the last event for this frame, if any.
    if let Some(event) = validation_events.read().last() {
        info!("[Validation] Triggered for SMILES: '{}'", event.0);
        let smiles = event.0.clone();

        // If a validation task is already running, cancel it.
        for (entity, mut task_component) in &mut existing_task_query {
            info!("[Validation] Cancelling previous validation task.");
            // Detaching the task handle prevents the result from being sent.
            if let Some(_old_task) = Option::take(&mut task_component.0) {
                // `old_task` is now owned by this scope. When the scope ends,
                // `old_task` is dropped, which signals bevy_tasks to cancel it.
                // We don't need to call any method on it.
            }
            commands.entity(entity).despawn();
            // commands.entity(entity).remove::<ValidationTask>();
        }

        // Use Bevy's async task pool to run the script in a background thread.
        let thread_pool = AsyncComputeTaskPool::get();
        let task = thread_pool.spawn(async move {
            if smiles.is_empty() {
                return SMILESValidationResult(Err("SMILES string is empty.".to_string()));
            }
            let output = std::process::Command::new("python")
                .arg("auto_typer.py")
                .arg("--smiles")
                .arg(&smiles)
                .output();

            info!("[Validation Task] Python script finished for '{}'.", smiles);

            match output {
                Ok(output) if output.status.success() => SMILESValidationResult(Ok(())),
                Ok(output) => {
                    let err_str = String::from_utf8_lossy(&output.stderr).to_string();
                    SMILESValidationResult(Err(err_str))
                }
                Err(e) => SMILESValidationResult(Err(e.to_string())),
            }
        });

        // Spawn a new entity to hold the task.
        commands.spawn(ValidationTask(Some(task)));
    }
}

fn handle_validation_result(
    mut commands: Commands,
    mut tasks: Query<(Entity, &mut ValidationTask)>,
    mut result_writer: EventWriter<SMILESValidationResult>,
) {
    for (entity, mut task_component) in &mut tasks {
        if let Some(task) = task_component.0.as_mut() {
            if let Some(result) = future::block_on(future::poll_once(task)) {
                info!(
                    "[Validation] Task complete. Sending result event: {:?}",
                    result
                );
                result_writer.write(result);
                commands.entity(entity).despawn();
            }
        }
    }
}

// System 1: Listens for external requests and runs the Python script.
fn trigger_molecule_generation_on_enter(
    mut spawn_events: EventReader<SpawnMoleculeFromSMILESEvent>,
    mut writer: EventWriter<SpawnMoleculeFromJsonEvent>,
) {
    for event in spawn_events.read() {
        let smiles = &event.0;
        info!("Received request to spawn molecule from SMILES: {}", smiles);
        let output = match std::process::Command::new("python")
            .arg("auto_typer.py")
            .arg("--smiles")
            .arg(smiles)
            .output()
        {
            Ok(out) => out,
            Err(e) => {
                error!("Failed to execute python script: {}", e);
                continue;
            }
        };
        if output.status.success() {
            writer.write(SpawnMoleculeFromJsonEvent(
                String::from_utf8_lossy(&output.stdout).to_string(),
            ));
        } else {
            error!(
                "Molecule generation script failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
    }
}

// System 2: Clears the board for the new molecule.
fn despawn_previous_molecule(
    mut commands: Commands,
    mut connectivity: ResMut<SystemConnectivity>,
    mut atom_id_map: ResMut<AtomIdMap>,
    // Also despawn any existing bond visuals
    bond_vis_query: Query<Entity, With<BondVisualization>>,
    molecule_container_query: Query<Entity, With<MoleculeContainer>>,
) {
    for entity in &molecule_container_query {
        commands.entity(entity).despawn();
    }
    for entity in &bond_vis_query {
        commands.entity(entity).despawn();
    }
    *connectivity = SystemConnectivity::default();
    *atom_id_map = AtomIdMap::default();
}

// System 3: Spawns the new molecule from the generated JSON.
fn spawn_new_molecule(
    mut events: EventReader<SpawnMoleculeFromJsonEvent>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut connectivity: ResMut<SystemConnectivity>,
    mut atom_id_map: ResMut<AtomIdMap>,
    shared_handles: Res<SharedAssetHandles>,
    force_field: Res<ForceField>,
    atom_query: Query<&Atom>,
) {
    for event in events.read() {
        let config: MoleculeConfig = match serde_json::from_str(&event.0) {
            Ok(c) => c,
            Err(e) => {
                error!("Failed to parse JSON from script: {}", e);
                continue;
            }
        };

        info!("Spawning new molecule: {}", config.name);
        let mut id_to_entity_map = HashMap::new();
        // The container needs the full set of visibility components to be a valid parent.
        let molecule_entity = commands
            .spawn((
                MoleculeContainer,
                Name::new(config.name.clone()),
                Transform::default(),
                GlobalTransform::default(),
                Visibility::default(),
                InheritedVisibility::default(),
                ViewVisibility::default(),
            ))
            .id();

        commands.entity(molecule_entity).with_children(|parent| {
            for atom_spec in &config.atoms {
                let radius = match atom_spec.element.as_str() {
                    "C" => 0.06,
                    "O" => 0.05,
                    "H" => 0.03,
                    "N" => 0.055,
                    _ => 0.045,
                };
                let color = match atom_spec.element.as_str() {
                    "C" => Color::srgb(0.2, 0.2, 0.2),
                    "O" => Color::srgb(1.0, 0.1, 0.1),
                    "H" => Color::srgb(0.9, 0.9, 0.9),
                    "N" => Color::srgb(0.1, 0.1, 1.0),
                    _ => Color::srgb(1.0, 0.2, 0.8),
                };

                // Individual atoms do not need Inherited/View visibility if they are children
                let entity = parent
                    .spawn((
                        Atom {
                            type_name: atom_spec.type_name.clone(),
                        },
                        Force::default(),
                        Velocity(Vec3::ZERO),
                        Acceleration(Vec3::ZERO),
                        Mesh3d(meshes.add(Sphere::new(radius))),
                        MeshMaterial3d(materials.add(color)),
                        Transform::from_translation(Vec3::from(atom_spec.pos)),
                        GlobalTransform::default(),
                    ))
                    .id();
                id_to_entity_map.insert(atom_spec.id.clone(), entity);
                atom_id_map
                    .entity_to_id
                    .insert(entity, atom_spec.id.clone());
            }
        });

        for bond_spec in &config.bonds {
            let entity1 = id_to_entity_map[&bond_spec.atoms[0]];
            let entity2 = id_to_entity_map[&bond_spec.atoms[1]];
            let order = match bond_spec.order.as_str() {
                "Double" => BondOrder::Double,
                "Triple" => BondOrder::Triple, // You need to add this to the BondOrder enum
                _ => BondOrder::Single,
            };
            connectivity.bonds.push(Bond {
                a: entity1,
                b: entity2,
                order,
            });
        }

        // Spawn bond visuals
        for bond in &connectivity.bonds {
            let Ok([atom1, atom2]) = atom_query.get_many([bond.a, bond.b]) else {
                continue;
            };
            let key = (atom1.type_name.clone(), atom2.type_name.clone(), bond.order);

            let material_handle = if force_field.bond_params.contains_key(&key) {
                shared_handles.bond_material.clone()
            } else {
                warn!(
                    "Found undefined bond between types: {} and {}",
                    atom1.type_name, atom2.type_name
                );
                shared_handles.undefined_bond_material.clone()
            };
            let num_strands = match bond.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
            };
            for i in 0..num_strands {
                commands.spawn((
                    Mesh3d(shared_handles.bond_mesh.clone()),
                    // Use the material we just selected
                    MeshMaterial3d(material_handle.clone()),
                    Transform::default(),
                    GlobalTransform::default(),
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
}

// System 4: Builds angles, dihedrals, and exclusions after spawning.
fn build_derived_connectivity(
    mut commands: Commands,
    mut connectivity: ResMut<SystemConnectivity>,
    // This event from interaction.rs signals other systems that a full rebuild happened.
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
) {
    // Only run if there are bonds to process.
    if connectivity.bonds.is_empty() {
        return;
    }

    info!("Building derived connectivity for new molecule...");
    let mut adjacency: HashMap<Entity, Vec<Entity>> = HashMap::new();
    for bond in &connectivity.bonds {
        adjacency.entry(bond.a).or_default().push(bond.b);
        adjacency.entry(bond.b).or_default().push(bond.a);
    }

    // `clone()` is needed here to satisfy the borrow checker, as we modify `connectivity` inside the loop.
    for (center_atom, neighbors) in adjacency.clone() {
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
    connectivity.dihedrals.append(&mut new_dihedrals);

    let mut one_two = HashSet::new();
    for bond in &connectivity.bonds {
        one_two.insert(if bond.a < bond.b {
            (bond.a, bond.b)
        } else {
            (bond.b, bond.a)
        });
    }
    let mut one_three = HashSet::new();
    for angle in &connectivity.angles {
        one_three.insert(if angle.a < angle.b {
            (angle.a, angle.b)
        } else {
            (angle.b, angle.a)
        });
    }
    commands.insert_resource(ExcludedPairs { one_two, one_three });

    info!(
        "Connectivity generation complete. Angles: {}, Dihedrals: {}",
        connectivity.angles.len(),
        connectivity.dihedrals.len()
    );
    rebuild_writer.write(RebuildConnectivityEvent);
}
