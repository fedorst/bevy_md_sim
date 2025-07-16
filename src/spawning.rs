// src/spawning.rs

use crate::components::*;
use crate::components::{Molecule, Solvent};
use crate::config::MoleculeConfig;
use crate::interaction::RebuildConnectivityEvent;
use crate::resources::{AtomIdMap, Bond, BondOrder, SimulationBox, SystemConnectivity};
use crate::spawning_utils::get_atom_visuals;
use bevy::prelude::*;
use bevy::tasks::{AsyncComputeTaskPool, Task};
use futures_lite::future;
use rand::Rng;
use std::collections::{HashMap, HashSet};

#[derive(Event, Debug)]
pub struct ClearAllMoleculesEvent;

/// An event sent to request a new molecule from a SMILES string.
#[derive(Event, Debug)]
pub struct SpawnMoleculeFromSMILESEvent(pub String, pub Option<IVec3>);

#[derive(Event, Debug)]
pub struct SpawnSolventEvent;

#[derive(Event, Debug)]
pub struct DespawnSolventEvent;

#[derive(Event, Debug)]
pub struct ValidateSMILESEvent(pub String);
/// Event sent by the async task with the validation result.
#[derive(Event, Debug)]
pub struct SMILESValidationResult(pub Result<(), String>);

/// An event sent when the python script succeeds, containing the molecule JSON.
#[derive(Event, Debug)]
pub struct SpawnMoleculeFromJsonEvent(pub String, pub Option<IVec3>);

#[derive(Component)]
struct ValidationTask(Option<Task<SMILESValidationResult>>);

pub struct SpawningPlugin;

impl Plugin for SpawningPlugin {
    fn build(&self, app: &mut App) {
        app.add_event::<ClearAllMoleculesEvent>()
            .add_event::<SpawnSolventEvent>()
            .add_event::<DespawnSolventEvent>()
            .add_event::<SpawnMoleculeFromSMILESEvent>()
            .add_event::<ValidateSMILESEvent>()
            .add_event::<SMILESValidationResult>()
            .add_event::<SpawnMoleculeFromJsonEvent>()
            .add_systems(
                Update,
                (
                    trigger_molecule_generation,
                    trigger_smiles_validation,
                    handle_validation_result,
                    clear_all_molecules_on_event,
                    spawn_molecules_from_json.run_if(on_event::<SpawnMoleculeFromJsonEvent>),
                    handle_spawn_solvent_event.run_if(on_event::<SpawnSolventEvent>),
                    handle_despawn_solvent_event.run_if(on_event::<DespawnSolventEvent>),
                ),
            );
        #[cfg(target_arch = "wasm32")]
        app.add_systems(Update, handle_molecule_generation_task);
    }
}

fn clear_all_molecules_on_event(
    mut commands: Commands,
    mut events: EventReader<ClearAllMoleculesEvent>,
    mut connectivity: ResMut<SystemConnectivity>,
    mut atom_id_map: ResMut<AtomIdMap>,
    bond_vis_query: Query<Entity, With<BondVisualization>>,
    molecule_query: Query<Entity, With<Molecule>>,
) {
    if events.read().last().is_none() {
        return;
    }
    info!("Clearing all molecules from the simulation.");
    for entity in &molecule_query {
        commands.entity(entity).despawn();
    }
    for entity in &bond_vis_query {
        commands.entity(entity).despawn();
    }
    *connectivity = SystemConnectivity::default();
    *atom_id_map = AtomIdMap::default();
}

fn handle_spawn_solvent_event(
    mut commands: Commands,
    mut events: EventReader<SpawnSolventEvent>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut connectivity: ResMut<SystemConnectivity>,
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
    sim_box_opt: Option<Res<SimulationBox>>,
    existing_atoms: Query<&Transform, With<Atom>>,
) {
    if events.read().last().is_none() {
        return;
    }
    let Some(sim_box) = sim_box_opt else {
        warn!("Cannot spawn solvent, no SimulationBox resource found.");
        return;
    };
    let output = match std::process::Command::new("python")
        .arg("src_py/auto_typer.py")
        .arg("--smiles")
        .arg("O")
        .output()
    {
        Ok(out) if out.status.success() => out,
        _ => {
            error!("Solvent generation script failed");
            return;
        }
    };
    let Ok(config) = serde_json::from_slice::<MoleculeConfig>(&output.stdout) else {
        return;
    };
    info!("Spawning solvent...");
    let grid_density = 5;
    let spacing = sim_box.size / grid_density as f32;
    let grid_offset = -(sim_box.size / 2.0) + (spacing / 2.0);
    let overlap_cutoff_sq = 0.3 * 0.3;

    let mut all_atom_positions: Vec<Vec3> = existing_atoms.iter().map(|t| t.translation).collect();

    for z in 0..grid_density {
        for y in 0..grid_density {
            for x in 0..grid_density {
                let molecule_offset =
                    Vec3::new(x as f32, y as f32, z as f32) * spacing + grid_offset;
                let is_overlapping = all_atom_positions
                    .iter()
                    .any(|pos| molecule_offset.distance_squared(*pos) < overlap_cutoff_sq);

                if !is_overlapping {
                    let mut id_to_entity_map = HashMap::new();
                    let mut spawned_atom_positions = Vec::new();
                    commands
                        .spawn((
                            Molecule,
                            Solvent,
                            (
                                Transform::default(),
                                GlobalTransform::default(),
                                Visibility::default(),
                                InheritedVisibility::default(),
                                ViewVisibility::default(),
                            ),
                        ))
                        .with_children(|parent| {
                            for atom_spec in &config.atoms {
                                let pos = Vec3::from(atom_spec.pos) + molecule_offset;
                                spawned_atom_positions.push(pos);
                                let (radius, color) = get_atom_visuals(&atom_spec.element);
                                let entity = parent
                                    .spawn((
                                        Atom {
                                            type_name: atom_spec.type_name.clone(),
                                        },
                                        Solvent,
                                        Force::default(),
                                        Velocity(Vec3::ZERO),
                                        Acceleration(Vec3::ZERO),
                                        Mesh3d(meshes.add(Sphere::new(radius))),
                                        MeshMaterial3d(materials.add(color)),
                                        Transform::from_translation(pos),
                                    ))
                                    .id();
                                id_to_entity_map.insert(atom_spec.id.clone(), entity);
                            }
                        });
                    all_atom_positions.extend(spawned_atom_positions);
                    for bond_spec in &config.bonds {
                        let entity1 = id_to_entity_map[&bond_spec.atoms[0]];
                        let entity2 = id_to_entity_map[&bond_spec.atoms[1]];
                        connectivity.bonds.push(Bond {
                            a: entity1,
                            b: entity2,
                            order: BondOrder::Single,
                        });
                    }
                }
            }
        }
    }
    rebuild_writer.write(RebuildConnectivityEvent);
}

fn handle_despawn_solvent_event(
    mut commands: Commands,
    mut events: EventReader<DespawnSolventEvent>,
    mut connectivity: ResMut<SystemConnectivity>,
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
    // Query for parent Molecule entities that are solvents
    solvent_molecule_query: Query<Entity, (With<Molecule>, With<Solvent>)>,
    // Query for child Atom entities that are solvents
    solvent_atom_query: Query<Entity, (With<Atom>, With<Solvent>)>,
) {
    if events.read().last().is_none() {
        return;
    }

    info!("Despawning solvent...");
    let solvent_atoms: HashSet<Entity> = solvent_atom_query.iter().collect();

    if solvent_atoms.is_empty() {
        return;
    }

    // 1. Remove all bonds connected to any solvent atom.
    connectivity
        .bonds
        .retain(|bond| !solvent_atoms.contains(&bond.a) && !solvent_atoms.contains(&bond.b));

    // 2. Despawn all solvent molecule entities, which will despawn their child atoms.
    for entity in &solvent_molecule_query {
        commands.entity(entity).despawn();
    }

    // 3. Trigger a full rebuild of visuals and derived connectivity.
    rebuild_writer.write(RebuildConnectivityEvent);
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
                .arg("src_py/auto_typer.py")
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
#[cfg(not(target_arch = "wasm32"))]
fn trigger_molecule_generation(
    mut spawn_events: EventReader<SpawnMoleculeFromSMILESEvent>,
    mut writer: EventWriter<SpawnMoleculeFromJsonEvent>,
) {
    for event in spawn_events.read() {
        let smiles = &event.0;
        let grid_info = event.1;
        info!("Received request to spawn molecule from SMILES: {}", smiles);
        let output = match std::process::Command::new("python")
            .arg("src_py/auto_typer.py")
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
                grid_info,
            ));
        } else {
            error!(
                "Molecule generation script failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
    }
}

#[cfg(target_arch = "wasm32")]
fn trigger_molecule_generation(
    mut commands: Commands,
    mut spawn_events: ResMut<Events<SpawnMoleculeFromSMILESEvent>>,
) {
    for event in spawn_events.drain() {
        let smiles = event.0.clone();
        let grid_info = event.1;

        let task = AsyncComputeTaskPool::get().spawn(async move {
            let api_url = "https://md.fedor.ee/api/generate_molecule";

            let client = reqwest::Client::new();
            let response = client
                .post(api_url)
                .json(&serde_json::json!({ "smiles": smiles }))
                .send()
                .await;

            match response {
                Ok(res) if res.status().is_success() => match res.text().await {
                    Ok(text) => Some(SpawnMoleculeFromJsonEvent(text, grid_info)),
                    Err(_) => None,
                },
                _ => None,
            }
        });

        // We need a component to hold the task and a system to poll it
        commands.spawn(MoleculeGenerationTask(task));
    }
}

#[cfg(target_arch = "wasm32")]
#[derive(Component)]
struct MoleculeGenerationTask(Task<Option<SpawnMoleculeFromJsonEvent>>);

#[cfg(target_arch = "wasm32")]
fn handle_molecule_generation_task(
    mut commands: Commands,
    mut tasks: Query<(Entity, &mut MoleculeGenerationTask)>,
    mut writer: EventWriter<SpawnMoleculeFromJsonEvent>,
) {
    for (entity, mut task) in &mut tasks {
        if let Some(Some(event)) =
            futures_lite::future::block_on(futures_lite::future::poll_once(&mut task.0))
        {
            info!(">>> [ASYNC_HANDLER] Task finished. Firing SpawnMoleculeFromJsonEvent.");
            writer.write(event);
            commands.entity(entity).despawn();
        }
    }
}

// System 3: Spawns the new molecule from the generated JSON.
fn spawn_molecules_from_json(
    mut events: EventReader<SpawnMoleculeFromJsonEvent>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut connectivity: ResMut<SystemConnectivity>,
    mut atom_id_map: ResMut<AtomIdMap>,
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
    existing_atoms: Query<&Transform, With<Atom>>,
) {
    for event in events.read() {
        info!(">>> [SPAWNER] Saw SpawnMoleculeFromJsonEvent. Spawning molecule NOW.");
        let config: MoleculeConfig = match serde_json::from_str(&event.0) {
            Ok(c) => c,
            Err(e) => {
                error!("Failed to parse JSON from script: {}", e);
                continue;
            }
        };

        let mut offset = Vec3::ZERO;
        if !existing_atoms.is_empty() {
            // 1. Find the center of the new molecule we are about to spawn
            let mut new_molecule_center = Vec3::ZERO;
            for atom_spec in &config.atoms {
                new_molecule_center += Vec3::from(atom_spec.pos);
            }
            new_molecule_center /= config.atoms.len() as f32;

            // 2. Find the bounding box of all existing atoms
            let mut min_bound = Vec3::splat(f32::MAX);
            let mut max_bound = Vec3::splat(f32::MIN);
            for transform in &existing_atoms {
                min_bound = min_bound.min(transform.translation);
                max_bound = max_bound.max(transform.translation);
            }
            let existing_size = (max_bound - min_bound).length();

            // 3. Create an offset vector along a random axis
            let mut rng = rand::rng();
            let mut random_axis = Vec3::ZERO;
            match rng.random_range(0..6) {
                0 => random_axis.x = 1.0,
                1 => random_axis.y = 1.0,
                2 => random_axis.z = 1.0,
                3 => random_axis.x = -1.0,
                4 => random_axis.y = -1.0,
                _ => random_axis.z = -1.0,
            };

            // The offset is the size of the existing system, plus a buffer,
            // translated away from the center of the new molecule.
            offset = (random_axis * (existing_size + 0.3)) - new_molecule_center;
        }

        let grid_size = event.1.unwrap_or(IVec3::ONE);
        let spacing = if event.1.is_some() { 1.2 } else { 2.4 };
        let box_size = if event.1.is_some() {
            grid_size.as_vec3() * spacing
        } else {
            Vec3::splat(spacing)
        };
        commands.insert_resource(SimulationBox { size: box_size });

        info!("Spawning {} molecule(s)...", config.name);
        let mut id_to_entity_map = HashMap::new();
        let molecule_entity = commands
            .spawn((
                Molecule,
                Name::new(config.name.clone()),
                (
                    Transform::default(),
                    GlobalTransform::default(),
                    Visibility::default(),
                    InheritedVisibility::default(),
                    ViewVisibility::default(),
                ),
            ))
            .id();
        commands.entity(molecule_entity).with_children(|parent| {
            for atom_spec in &config.atoms {
                let (radius, color) = get_atom_visuals(&atom_spec.element);
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
                        Transform::from_translation(Vec3::from(atom_spec.pos) + offset),
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
                "Triple" => BondOrder::Triple,
                _ => BondOrder::Single,
            };
            connectivity.bonds.push(Bond {
                a: entity1,
                b: entity2,
                order,
            });
        }
    }
    rebuild_writer.write(RebuildConnectivityEvent);
}
