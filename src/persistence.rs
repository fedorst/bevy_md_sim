// src/persistence.rs

use crate::components::{Acceleration, Atom, Force, Molecule, Velocity};
use crate::interaction::RebuildConnectivityEvent;
use crate::resources::{AtomIdMap, Bond, BondOrder, LastSaveTime, SystemConnectivity};
#[cfg(target_arch = "wasm32")]
use crate::spawning::SpawnMoleculeFromJsonEvent;
use crate::spawning_utils::get_atom_visuals;
use bevy::prelude::*; // Make sure Time is imported
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::{Duration, SystemTime};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
extern "C" {
    fn js_save_to_localstorage(key: &str, value: &str);
    fn js_load_from_localstorage(key: &str) -> Option<String>;
    fn js_get_save_timestamp() -> Option<String>;
}

const SAVE_FILE_PATH: &str = "simulation_state.json";
const LOCALSTORAGE_KEY: &str = "md_simulation_save_state";

fn create_save_state(
    atom_query: &Query<(&Atom, &Transform, &Velocity, Entity)>,
    connectivity: &Res<SystemConnectivity>,
    atom_id_map: &Res<AtomIdMap>,
) -> SimulationSaveState {
    info!("Serializing simulation state...");

    let atom_states: Vec<AtomState> = atom_query
        .iter()
        .filter_map(|(atom, transform, velocity, entity)| {
            atom_id_map.entity_to_id.get(&entity).map(|id| {
                let element = atom.type_name.chars().next().unwrap_or('?').to_string();
                AtomState {
                    id: id.clone(),
                    type_name: atom.type_name.clone(),
                    element,
                    position: transform.translation.to_array(),
                    velocity: velocity.0.to_array(),
                }
            })
        })
        .collect();

    let saved_bonds: Vec<SavedBond> = connectivity
        .bonds
        .iter()
        .filter_map(|bond| {
            let id1 = atom_id_map.entity_to_id.get(&bond.a)?;
            let id2 = atom_id_map.entity_to_id.get(&bond.b)?;
            Some(SavedBond {
                atom_ids: [id1.clone(), id2.clone()],
                order: bond.order,
            })
        })
        .collect();

    SimulationSaveState {
        atoms: atom_states,
        bonds: saved_bonds,
    }
}

#[cfg(target_arch = "wasm32")]
fn save_state_to_localstorage_on_event(
    mut events: EventReader<SaveStateEvent>,
    atom_query: Query<(&Atom, &Transform, &Velocity, Entity)>,
    connectivity: Res<SystemConnectivity>,
    atom_id_map: Res<AtomIdMap>,
) {
    if events.read().last().is_none() {
        return;
    }

    let save_state = create_save_state(&atom_query, &connectivity, &atom_id_map);

    match serde_json::to_string(&save_state) {
        Ok(json_string) => {
            js_save_to_localstorage(LOCALSTORAGE_KEY, &json_string);
            info!("Successfully saved state to LocalStorage.");
        }
        Err(e) => error!("Failed to serialize state for LocalStorage: {}", e),
    }
}

#[cfg(target_arch = "wasm32")]
fn load_state_from_localstorage(
    mut commands: Commands,
    mut load_events: EventReader<LoadStateEvent>,
    mut connectivity: ResMut<SystemConnectivity>,
    mut atom_id_map: ResMut<AtomIdMap>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
    molecule_query: Query<Entity, With<Molecule>>,
) {
    if load_events.read().last().is_none() {
        return;
    }

    info!("[LOAD] Attempting to load state from LocalStorage.");
    if let Some(json_string) = js_load_from_localstorage(LOCALSTORAGE_KEY) {
        let Ok(save_state) = serde_json::from_str::<SimulationSaveState>(&json_string) else {
            error!("[LOAD] Failed to deserialize save state from LocalStorage.");
            return;
        };

        // 1. Clear the stage.
        for entity in &molecule_query {
            commands.entity(entity).despawn();
        }
        *connectivity = SystemConnectivity::default();
        *atom_id_map = AtomIdMap::default();

        // 2. Spawn atoms from the save state using the correct 0.16 component tuple syntax.
        let mut id_to_entity_map = HashMap::new();

        commands
            .spawn((
                Molecule,
                Name::new("Loaded Molecule"),
                (
                    Transform::default(),
                    GlobalTransform::default(),
                    Visibility::default(),
                    InheritedVisibility::default(),
                    ViewVisibility::default(),
                ),
            ))
            .with_children(|parent| {
                for atom_state in &save_state.atoms {
                    let (radius, color) = get_atom_visuals(&atom_state.element);
                    let entity = parent
                        .spawn((
                            Atom {
                                type_name: atom_state.type_name.clone(),
                            },
                            Force::default(),
                            Velocity(Vec3::from(atom_state.position)),
                            Acceleration(Vec3::ZERO),
                            Mesh3d(meshes.add(Sphere::new(radius))),
                            MeshMaterial3d(materials.add(color)),
                            Transform::from_translation(Vec3::from(atom_state.position)),
                        ))
                        .id();
                    id_to_entity_map.insert(atom_state.id.clone(), entity);
                    atom_id_map
                        .entity_to_id
                        .insert(entity, atom_state.id.clone());
                }
            });

        // 3. Rebuild bond connectivity
        for saved_bond in &save_state.bonds {
            let (Some(entity1), Some(entity2)) = (
                id_to_entity_map.get(&saved_bond.atom_ids[0]),
                id_to_entity_map.get(&saved_bond.atom_ids[1]),
            ) else {
                continue;
            };
            connectivity.bonds.push(Bond {
                a: *entity1,
                b: *entity2,
                order: saved_bond.order,
            });
        }

        // 4. Trigger a full rebuild of visuals
        rebuild_writer.write(RebuildConnectivityEvent);
        info!("[LOAD] Successfully loaded and spawned state from LocalStorage.");
    } else {
        warn!("[LOAD] No save data found in LocalStorage.");
    }
}

#[cfg(target_arch = "wasm32")]
fn update_save_time_display_wasm(
    mut last_save_time: ResMut<LastSaveTime>,
    time: Res<Time>,
    mut timer: Local<Timer>,
) {
    timer.set_duration(Duration::from_secs(2)); // Check every 2 seconds
    timer.tick(time.delta());

    if !timer.finished() {
        return;
    }

    if let Some(timestamp_str) = js_get_save_timestamp() {
        // Here we would parse the ISO string and format it nicely.
        // For now, a simple display is fine.
        // A proper implementation would use a chrono-like library.
        last_save_time.display_text = format!("Last save: {}", &timestamp_str[11..19]);
    } else {
        last_save_time.display_text = "No save file".to_string();
    }
}

pub struct PersistencePlugin;

impl Plugin for PersistencePlugin {
    fn build(&self, app: &mut App) {
        app.add_event::<SaveStateEvent>()
            .add_event::<LoadStateEvent>();

        #[cfg(not(target_arch = "wasm32"))]
        app.add_systems(
            Update,
            (
                save_state_to_file_on_event,
                load_state_from_file_on_event,
                update_save_time_display,
            ),
        );

        // Systems for web LocalStorage-based persistence
        #[cfg(target_arch = "wasm32")]
        app.add_systems(
            Update,
            (
                save_state_to_localstorage_on_event,
                load_state_from_localstorage,
                update_save_time_display_wasm,
            ),
        );
    }
}

// --- Events ---
#[derive(Event)]
pub struct SaveStateEvent;

#[derive(Event)]
pub struct LoadStateEvent;

// --- Data Structures for Serialization
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct AtomState {
    pub id: String,
    pub type_name: String,
    pub element: String,
    pub position: [f32; 3],
    pub velocity: [f32; 3],
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct SimulationSaveState {
    pub atoms: Vec<AtomState>,
    pub bonds: Vec<SavedBond>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SavedBond {
    pub atom_ids: [String; 2],
    pub order: BondOrder,
}

// --- Systems ---

#[cfg(not(target_arch = "wasm32"))]
fn save_state_to_file_on_event(
    mut events: EventReader<SaveStateEvent>,
    atom_query: Query<(&Atom, &Transform, &Velocity, Entity)>,
    connectivity: Res<SystemConnectivity>,
    atom_id_map: Res<AtomIdMap>,
) {
    if events.read().last().is_none() {
        return;
    }

    let save_state = create_save_state(&atom_query, &connectivity, &atom_id_map);

    match serde_json::to_string_pretty(&save_state) {
        Ok(json_string) => {
            if let Err(e) = std::fs::write(SAVE_FILE_PATH, json_string) {
                error!("[SAVE] FAILED to write save file: {}", e);
            } else {
                info!("[SAVE] Successfully wrote state to {}.", SAVE_FILE_PATH);
            }
        }
        Err(e) => error!("[SAVE] FAILED to serialize save state: {}", e),
    }
}

fn load_state_from_file_on_event(
    mut events: EventReader<LoadStateEvent>,
    mut commands: Commands,
    mut connectivity: ResMut<SystemConnectivity>,
    mut atom_id_map: ResMut<AtomIdMap>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
    molecule_query: Query<Entity, With<Molecule>>,
) {
    if let Some(_event) = events.read().last() {
        let Ok(json_string) = std::fs::read_to_string(SAVE_FILE_PATH) else {
            error!("[LOAD] FAILED to read save file: {}", SAVE_FILE_PATH);
            return;
        };
        let Ok(save_state) = serde_json::from_str::<SimulationSaveState>(&json_string) else {
            error!(
                "[LOAD] FAILED to deserialize save state from {}",
                SAVE_FILE_PATH
            );
            return;
        };
        info!(
            "[LOAD] Successfully deserialized save state with {} atoms and {} bonds.",
            save_state.atoms.len(),
            save_state.bonds.len()
        );

        // 1. Clear the stage.
        info!("[LOAD] Despawning all existing molecules...");
        for entity in &molecule_query {
            commands.entity(entity).despawn();
        }
        *connectivity = SystemConnectivity::default();
        *atom_id_map = AtomIdMap::default();

        // 2. Spawn atoms from the save file.
        let mut id_to_entity_map = HashMap::new();
        let molecule_entity = commands
            .spawn((
                Molecule,
                Name::new("Loaded State"),
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
            for atom_state in &save_state.atoms {
                let (radius, color) = get_atom_visuals(&atom_state.element);
                let entity = parent
                    .spawn((
                        Atom {
                            type_name: atom_state.type_name.clone(),
                        },
                        Force::default(),
                        Velocity(Vec3::from(atom_state.velocity)),
                        Acceleration(Vec3::ZERO),
                        Mesh3d(meshes.add(Sphere::new(radius))),
                        MeshMaterial3d(materials.add(color)),
                        Transform::from_translation(Vec3::from(atom_state.position)),
                    ))
                    .id();
                id_to_entity_map.insert(atom_state.id.clone(), entity);
                atom_id_map
                    .entity_to_id
                    .insert(entity, atom_state.id.clone());
            }
        });
        info!("[LOAD] Spawned {} new atoms.", id_to_entity_map.len());

        // 3. Rebuild bond connectivity.
        for saved_bond in &save_state.bonds {
            let Some(entity1) = id_to_entity_map.get(&saved_bond.atom_ids[0]) else {
                continue;
            };
            let Some(entity2) = id_to_entity_map.get(&saved_bond.atom_ids[1]) else {
                continue;
            };
            connectivity.bonds.push(Bond {
                a: *entity1,
                b: *entity2,
                order: saved_bond.order,
            });
        }
        info!(
            "[LOAD] Recreated {} bonds in connectivity resource.",
            connectivity.bonds.len()
        );

        // 4. Trigger a full rebuild of visuals and derived connectivity.
        rebuild_writer.write(RebuildConnectivityEvent);
        info!("[LOAD] Finished loading state. Firing RebuildConnectivityEvent.");
    }
}

fn update_save_time_display(
    mut last_save_time: ResMut<LastSaveTime>,
    // This system runs on a timer to avoid checking the file system every single frame.
    time: Res<Time>,
    mut timer: Local<Timer>,
) {
    timer.set_duration(Duration::from_secs(1)); // Check once per second
    timer.tick(time.delta());
    if !timer.finished() {
        return;
    }

    if let Ok(metadata) = std::fs::metadata(SAVE_FILE_PATH) {
        if let Ok(modified_time) = metadata.modified() {
            if let Ok(elapsed) = SystemTime::now().duration_since(modified_time) {
                // THE FIX: Format the duration into a human-readable string.
                last_save_time.display_text = format_time_ago(elapsed.as_secs());
            } else {
                last_save_time.display_text = "Last save: now".to_string();
            }
        } else {
            last_save_time.display_text = "Save time invalid".to_string();
        }
    } else {
        last_save_time.display_text = "No save file".to_string();
    }
}

// THE FIX: A new helper function to format seconds into a nice string.
fn format_time_ago(seconds: u64) -> String {
    if seconds < 2 {
        return "Last save: just now".to_string();
    }
    if seconds < 60 {
        return format!("Last save: {}s ago", seconds);
    }

    let minutes = seconds / 60;
    if minutes < 60 {
        return format!("Last save: {}m ago", minutes);
    }

    let hours = minutes / 60;
    if hours < 24 {
        return format!("Last save: {}h ago", hours);
    }

    let days = hours / 24;
    format!("Last save: {}d ago", days)
}
