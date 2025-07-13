use bevy::prelude::*;
use serde::Deserialize;
use std::collections::VecDeque;
use std::collections::{HashMap, HashSet};
use std::time::Duration;

#[derive(Resource, Default)]
pub struct LastClick {
    pub time: Option<Duration>,
    pub target: Option<Entity>,
}

#[derive(Resource)]
pub struct SimulationBox {
    pub size: Vec3,
}

#[derive(Resource)]
pub struct EnergyHistory {
    pub potential: VecDeque<(f64, f64)>, // (step, value)
    pub kinetic: VecDeque<(f64, f64)>,
    pub total: VecDeque<(f64, f64)>,
    pub capacity: usize,
}

impl Default for EnergyHistory {
    fn default() -> Self {
        let capacity = 500; // Store the last 500 steps
        Self {
            potential: VecDeque::with_capacity(capacity),
            kinetic: VecDeque::with_capacity(capacity),
            total: VecDeque::with_capacity(capacity),
            capacity,
        }
    }
}

// data structures used by resources
#[derive(Debug, Copy, Clone)]
pub struct Bond {
    pub a: Entity,
    pub b: Entity,
    pub order: BondOrder,
}

#[derive(Debug, Copy, Clone)]
pub struct Angle {
    pub a: Entity,
    pub center: Entity,
    pub b: Entity,
}
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

// resources

#[derive(Resource)]
pub struct SimulationParameters {
    pub dt: f32,
}

#[derive(Resource, Default)]
pub struct PauseMenuState {
    pub visible: bool,
}

#[derive(Resource, Deref, DerefMut)]
pub struct ForceMultiplier(pub f32);

#[derive(Resource)]
pub struct SharedAssetHandles {
    pub bond_mesh: Handle<Mesh>,
    pub bond_material: Handle<StandardMaterial>,
    pub undefined_bond_material: Handle<StandardMaterial>,
}

#[derive(Resource, Default, Deref, DerefMut)]
pub struct StepCount(pub u32);

#[derive(Resource, Default)]
pub struct DragState {
    pub dragged_entities: Vec<Entity>,
    pub initial_positions: Vec<Vec3>,
    pub initial_centroid: Option<Vec3>,
    pub plane: Option<InfinitePlane3d>,
    pub target_centroid: Option<Vec3>,
}

#[derive(Resource, Default, Debug, Clone, Copy)]
pub struct SystemEnergy {
    pub potential: f32,
    pub kinetic: f32,
    pub total: f32,
}

// --- Structs for Deserializing the Force Field File ---
#[derive(Deserialize, Debug, Clone)]
pub struct AtomTypeParam {
    pub element: String,
    pub mass: f32,
    pub charge: f32,
    pub sigma: f32,
    pub epsilon: f32,
}

#[derive(Deserialize, Debug)]
struct BondParam {
    types: [String; 2],
    order: String, // "Single" or "Double"
    k: f32,
    r0: f32,
}

#[derive(Deserialize, Debug)]
struct AngleParam {
    types: [String; 3],
    k: f32,
    theta0_deg: f32,
}

#[derive(Deserialize, Debug)]
struct DihedralParam {
    types: [String; 4],
    k: f32,
    n: i32,
    phi0_deg: f32,
}

#[derive(Deserialize, Debug)]
struct ForceFieldFile {
    atom_types: HashMap<String, AtomTypeParam>,
    bonds: Vec<BondParam>,
    angles: Vec<AngleParam>,
    dihedrals: Vec<DihedralParam>,
}

pub struct Dihedral {
    pub a: Entity,
    pub b: Entity,
    pub c: Entity,
    pub d: Entity,
}

#[derive(Resource, Default)]
pub struct SystemConnectivity {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>,
}

#[derive(Resource, Debug, Default)]
pub struct ForceField {
    pub atom_types: HashMap<String, AtomTypeParam>,
    pub bond_params: HashMap<(String, String, BondOrder), (f32, f32)>,
    pub angle_params: HashMap<(String, String, String), (f32, f32)>,
    pub dihedral_params: HashMap<(String, String, String, String), (f32, i32, f32)>,
    pub coulomb_k: f32,
}

impl ForceField {
    pub fn from_file(path: &str) -> Self {
        let file_contents = std::fs::read_to_string(path)
            .unwrap_or_else(|_| panic!("Failed to read force field file at {}", path));
        let ff_file: ForceFieldFile =
            serde_json::from_str(&file_contents).expect("Failed to parse force field JSON");
        let mut bond_params = HashMap::new();
        for p in ff_file.bonds {
            let order = match p.order.as_str() {
                "Double" => BondOrder::Double,
                "Triple" => BondOrder::Triple,
                _ => BondOrder::Single,
            };
            bond_params.insert((p.types[0].clone(), p.types[1].clone(), order), (p.k, p.r0));
            bond_params.insert((p.types[1].clone(), p.types[0].clone(), order), (p.k, p.r0));
        }

        let mut angle_params = HashMap::new();

        for p in ff_file.angles {
            let theta_rad = p.theta0_deg.to_radians();
            angle_params.insert(
                (p.types[0].clone(), p.types[1].clone(), p.types[2].clone()),
                (p.k, theta_rad),
            );
            angle_params.insert(
                (p.types[2].clone(), p.types[1].clone(), p.types[0].clone()),
                (p.k, theta_rad),
            );
        }

        let mut dihedral_params = HashMap::new();
        for p in ff_file.dihedrals {
            let phi0_rad = p.phi0_deg.to_radians();
            dihedral_params.insert(
                (
                    p.types[0].clone(),
                    p.types[1].clone(),
                    p.types[2].clone(),
                    p.types[3].clone(),
                ),
                (p.k, p.n, phi0_rad),
            );
            dihedral_params.insert(
                (
                    p.types[3].clone(),
                    p.types[2].clone(),
                    p.types[1].clone(),
                    p.types[0].clone(),
                ),
                (p.k, p.n, phi0_rad),
            );
        }

        Self {
            atom_types: ff_file.atom_types,
            bond_params,
            angle_params,
            dihedral_params,
            coulomb_k: 138.935458,
        }
    }
}

#[derive(Resource, Default)]
pub struct ExcludedPairs {
    pub one_two: HashSet<(Entity, Entity)>,   // Direct bonds (1-2)
    pub one_three: HashSet<(Entity, Entity)>, // Angle partners (1-3)
}

#[derive(Resource, Default, Deref, DerefMut)]
pub struct ActiveWallTime(pub f32);

#[derive(Resource, Default)]
pub struct Thermostat {
    pub target_temperature: f32, // in K
    pub tau: f32,                // coupling time in ps
}

#[derive(Resource, Default, Deref, DerefMut)]
pub struct CurrentTemperature(pub f32);

#[derive(Resource, Default, Deref, DerefMut)]
pub struct ThermostatScale(pub f32);

#[derive(Resource, Default)]
pub struct AtomIdMap {
    pub entity_to_id: HashMap<Entity, String>,
}
