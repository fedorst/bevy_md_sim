use bevy::prelude::*;
use std::collections::{HashMap, HashSet};

#[derive(Component, Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum AtomType {
    Oxygen,
    Hydrogen,
    Carbon,
}

impl AtomType {
    pub fn mass(&self) -> f32 {
        match self {
            AtomType::Oxygen => 15.999,
            AtomType::Hydrogen => 1.0,
            AtomType::Carbon => 12.011,
        }
    }
    pub fn charge(&self) -> f32 {
        match self {
            AtomType::Oxygen => -0.683,  // O in an alcohol
            AtomType::Hydrogen => 0.418, // H bonded to the O
            // We need to distinguish between carbons and their hydrogens.
            // For now, we'll use a simplified model and assign an "average" charge. TODO
            // A more advanced system would have more AtomTypes (e.g., C_CH3, H_CH3).
            AtomType::Carbon => -0.1,
        }
    }
    pub fn epsilon(&self) -> f32 {
        // OPLS-AA parameters (kJ/mol)
        match self {
            Self::Oxygen => 0.71128,
            Self::Hydrogen => 0.0, // Typically non-interacting for LJ
            Self::Carbon => 0.27614,
        }
    }
    pub fn sigma(&self) -> f32 {
        match self {
            Self::Oxygen => 0.3166,
            Self::Hydrogen => 0.0,
            Self::Carbon => 0.350,
        }
    }
}

#[derive(Component)]
pub struct BondVisualization {
    pub atom1: Entity,
    pub atom2: Entity,
}

#[derive(Resource, Default, PartialEq, Eq, Clone, Copy, Debug)]
pub struct StepSimulation(pub bool);

#[derive(Resource)]
pub struct SimulationParameters {
    pub dt: f32, // picoseconds
}

#[derive(Resource, Default)]
pub struct SimulationState {
    pub paused: bool,
}

#[derive(Resource, Default, Deref, DerefMut)]
pub struct StepCount(pub u32);

// A single bond between two atom entities.
#[derive(Debug, Copy, Clone)]
pub struct Bond {
    pub a: Entity,
    pub b: Entity,
}

// A single angle between three atom entities.
#[derive(Debug, Copy, Clone)]
pub struct Angle {
    pub a: Entity,
    pub center: Entity,
    pub b: Entity,
}

// A global resource that stores all bonds and angles in the system.
#[derive(Resource, Default)]
pub struct SystemConnectivity {
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
}

#[derive(Component, Deref, DerefMut, Default, Debug)]
pub struct Force(pub Vec3);

#[derive(Resource, Debug)]
pub struct ForceField {
    pub bond_params: HashMap<(AtomType, AtomType), (f32, f32)>, // (k, r0)
    pub angle_params: HashMap<(AtomType, AtomType, AtomType), (f32, f32)>, // (k, theta0)
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

        Self {
            bond_params,
            angle_params,
            coulomb_k: 138.935458,
        }
    }
}

#[derive(Component, Deref, DerefMut, Debug)]
pub struct Velocity(pub Vec3);

#[derive(Component, Deref, DerefMut)]
pub struct Acceleration(pub Vec3);

#[derive(Resource, Default)]
pub struct ExcludedPairs(pub HashSet<(Entity, Entity)>);

#[derive(Resource, Default, Deref, DerefMut)]
pub struct ActiveWallTime(pub f32);

pub struct CorePlugin;

impl Plugin for CorePlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<ForceField>()
            .init_resource::<SystemConnectivity>()
            .init_resource::<StepSimulation>()
            .insert_resource(SimulationParameters { dt: 1e-5 })
            .init_resource::<StepCount>()
            .init_resource::<SimulationState>()
            .init_resource::<ActiveWallTime>()
            .init_resource::<ExcludedPairs>();
    }
}
