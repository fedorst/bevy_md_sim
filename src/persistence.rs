// src/persistence.rs

use crate::components::{Acceleration, Atom, Force, Molecule, Velocity};
use crate::interaction::RebuildConnectivityEvent;
use crate::resources::{AtomIdMap, Bond, BondOrder, LastSaveTime, SystemConnectivity};
use bevy::prelude::*; // Make sure Time is imported
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::{Duration, SystemTime};

pub struct PersistencePlugin;
pub const SAVE_FILE_PATH: &str = "simulation_state.json";

impl Plugin for PersistencePlugin {
    fn build(&self, app: &mut App) {
        app.add_event::<SaveStateEvent>()
            .add_event::<LoadStateEvent>()
            .add_systems(
                Update,
                (
                    save_state_on_event,
                    load_state_on_event,
                    update_save_time_display,
                ),
            );
    }
}

// --- Events ---
#[derive(Event)]
pub struct SaveStateEvent(pub String); // Carries the filepath

#[derive(Event)]
pub struct LoadStateEvent(pub String); // Carries the filepath

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

fn save_state_on_event(
    mut events: EventReader<SaveStateEvent>,
    atom_query: Query<(&Atom, &Transform, &Velocity, Entity)>,
    connectivity: Res<SystemConnectivity>,
    atom_id_map: Res<AtomIdMap>,
) {
    if let Some(event) = events.read().last() {
        let filepath = &event.0;
        info!("[SAVE] Save event received for '{}'.", filepath);

        // 1. Convert atoms to serializable state.
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
        info!("[SAVE] Found and serialized {} atoms.", atom_states.len());

        // 2. Convert bonds to serializable state.
        let saved_bonds: Vec<SavedBond> = connectivity
            .bonds
            .iter()
            .filter_map(|bond| {
                let Some(id1) = atom_id_map.entity_to_id.get(&bond.a) else {
                    return None;
                };
                let Some(id2) = atom_id_map.entity_to_id.get(&bond.b) else {
                    return None;
                };
                Some(SavedBond {
                    atom_ids: [id1.clone(), id2.clone()],
                    order: bond.order,
                })
            })
            .collect();
        info!("[SAVE] Found and serialized {} bonds.", saved_bonds.len());

        let save_state = SimulationSaveState {
            atoms: atom_states,
            bonds: saved_bonds,
        };

        match serde_json::to_string_pretty(&save_state) {
            Ok(json_string) => {
                if let Err(e) = std::fs::write(filepath, json_string) {
                    error!("[SAVE] FAILED to write save file: {}", e);
                } else {
                    info!("[SAVE] Successfully wrote state to {}.", filepath);
                }
            }
            Err(e) => error!("[SAVE] FAILED to serialize save state: {}", e),
        }
    }
}

fn load_state_on_event(
    mut events: EventReader<LoadStateEvent>,
    mut commands: Commands,
    mut connectivity: ResMut<SystemConnectivity>,
    mut atom_id_map: ResMut<AtomIdMap>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut rebuild_writer: EventWriter<RebuildConnectivityEvent>,
    molecule_query: Query<Entity, With<Molecule>>,
) {
    if let Some(event) = events.read().last() {
        let filepath = &event.0;
        info!("[LOAD] Load event received for '{}'.", filepath);

        let Ok(json_string) = std::fs::read_to_string(filepath) else {
            error!("[LOAD] FAILED to read save file: {}", filepath);
            return;
        };
        let Ok(save_state) = serde_json::from_str::<SimulationSaveState>(&json_string) else {
            error!("[LOAD] FAILED to deserialize save state from {}", filepath);
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
                let radius = match atom_state.element.as_str() {
                    "C" => 0.06,
                    "O" => 0.05,
                    "H" => 0.03,
                    "N" => 0.055,
                    "S" => 0.1,
                    _ => 0.045,
                };
                let color = match atom_state.element.as_str() {
                    "C" => Color::srgb(0.2, 0.2, 0.2),
                    "O" => Color::srgb(1.0, 0.1, 0.1),
                    "H" => Color::srgb(0.9, 0.9, 0.9),
                    "N" => Color::srgb(0.1, 0.1, 1.0),
                    "S" => Color::srgb(1.0, 1.0, 0.0),
                    _ => Color::srgb(1.0, 0.2, 0.8),
                };

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
