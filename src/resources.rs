use crate::components::AtomType;
use bevy::prelude::*;
use std::collections::{HashMap, HashSet};

// data structures used by resources
#[derive(Debug, Copy, Clone)]
pub struct Bond {
    pub a: Entity,
    pub b: Entity,
}

#[derive(Debug, Copy, Clone)]
pub struct Angle {
    pub a: Entity,
    pub center: Entity,
    pub b: Entity,
}

// resources

#[derive(Resource, Default, PartialEq, Eq, Clone, Copy, Debug)]
pub struct StepSimulation(pub bool);

#[derive(Resource)]
pub struct SimulationParameters {
    pub dt: f32,
}

#[derive(Resource)]
pub struct MoleculeSelection(pub String);

#[derive(Resource, Default)]
pub struct SimulationState {
    pub paused: bool,
}

#[derive(Resource)]
pub struct SharedAssetHandles {
    pub bond_mesh: Handle<Mesh>,
    pub bond_material: Handle<StandardMaterial>,
}

#[derive(Resource, Default, Deref, DerefMut)]
pub struct StepCount(pub u32);

#[derive(Resource, Default, Debug, Clone, Copy)]
pub struct SystemEnergy {
    pub potential: f32,
    pub kinetic: f32,
    pub total: f32,
}

#[derive(Resource, Default)]
pub struct SystemConnectivity {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
}

#[derive(Resource, Debug)]
pub struct ForceField {
    pub bond_params: HashMap<(AtomType, AtomType), (f32, f32)>,
    pub angle_params: HashMap<(AtomType, AtomType, AtomType), (f32, f32)>,
    pub coulomb_k: f32,
}

impl Default for ForceField {
    fn default() -> Self {
        let mut bond_params = HashMap::new();
        // Insert params for both orders (C,O) and (O,C) for easier lookup
        // OPLS-AA parameters
        // C-C bond
        bond_params.insert((AtomType::Carbon, AtomType::Carbon), (224262.4, 0.1529));
        // C-H bond
        bond_params.insert((AtomType::Carbon, AtomType::Hydrogen), (284512.0, 0.1090));
        bond_params.insert((AtomType::Hydrogen, AtomType::Carbon), (284512.0, 0.1090));
        // C-O bond
        bond_params.insert((AtomType::Carbon, AtomType::Oxygen), (267776.0, 0.1410));
        bond_params.insert((AtomType::Oxygen, AtomType::Carbon), (267776.0, 0.1410));
        // O-H bond
        bond_params.insert((AtomType::Oxygen, AtomType::Hydrogen), (462750.4, 0.0945));
        bond_params.insert((AtomType::Hydrogen, AtomType::Oxygen), (462750.4, 0.0945));
        // C-N bond
        bond_params.insert((AtomType::Carbon, AtomType::Nitrogen), (334720.0, 0.147));
        bond_params.insert((AtomType::Nitrogen, AtomType::Carbon), (334720.0, 0.147));
        // N-H bond
        bond_params.insert((AtomType::Nitrogen, AtomType::Hydrogen), (389112.0, 0.101));
        bond_params.insert((AtomType::Hydrogen, AtomType::Nitrogen), (389112.0, 0.101));

        let mut angle_params = HashMap::new();
        // C-C-H angle
        angle_params.insert(
            (AtomType::Carbon, AtomType::Carbon, AtomType::Hydrogen),
            (292.88, 110.7f32.to_radians()),
        );
        angle_params.insert(
            (AtomType::Hydrogen, AtomType::Carbon, AtomType::Carbon),
            (292.88, 110.7f32.to_radians()),
        );
        // H-C-H angle
        angle_params.insert(
            (AtomType::Hydrogen, AtomType::Carbon, AtomType::Hydrogen),
            (276.144, 107.8f32.to_radians()),
        );
        // C-C-O angle
        angle_params.insert(
            (AtomType::Carbon, AtomType::Carbon, AtomType::Oxygen),
            (418.4, 108.5f32.to_radians()),
        );
        angle_params.insert(
            (AtomType::Oxygen, AtomType::Carbon, AtomType::Carbon),
            (418.4, 108.5f32.to_radians()),
        );
        // C-O-H angle
        angle_params.insert(
            (AtomType::Carbon, AtomType::Oxygen, AtomType::Hydrogen),
            (460.24, 108.5f32.to_radians()),
        );
        angle_params.insert(
            (AtomType::Hydrogen, AtomType::Oxygen, AtomType::Carbon),
            (460.24, 108.5f32.to_radians()),
        );

        // C-N-H angle
        angle_params.insert(
            (AtomType::Carbon, AtomType::Nitrogen, AtomType::Hydrogen),
            (292.88, 109.5f32.to_radians()),
        );
        angle_params.insert(
            (AtomType::Hydrogen, AtomType::Nitrogen, AtomType::Carbon),
            (292.88, 109.5f32.to_radians()),
        );
        // H-N-H angle
        angle_params.insert(
            (AtomType::Hydrogen, AtomType::Nitrogen, AtomType::Hydrogen),
            (292.88, 107.0f32.to_radians()),
        );
        // C-N-C angle
        angle_params.insert(
            (AtomType::Carbon, AtomType::Nitrogen, AtomType::Carbon),
            (418.4, 112.0f32.to_radians()),
        );

        Self {
            bond_params,
            angle_params,
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
