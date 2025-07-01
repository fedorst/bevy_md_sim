use bevy::prelude::*;

#[derive(Component, Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum AtomType {
    Oxygen,
    Hydrogen,
    Carbon,
    Nitrogen,
}

impl AtomType {
    pub fn mass(&self) -> f32 {
        match self {
            AtomType::Oxygen => 15.999,
            AtomType::Hydrogen => 1.0,
            AtomType::Carbon => 12.011,
            AtomType::Nitrogen => 14.007, // NEW
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
            // OPLS-AA charge for an amine nitrogen
            AtomType::Nitrogen => -1.05,
            // NOTE: The hydrogens on this nitrogen have a different charge
            // than the hydrogens on oxygen. This highlights the need for a
            // more advanced "atom type name" system later, but for now,
            // we will use the existing Hydrogen charge as an approximation.
        }
    }
    pub fn epsilon(&self) -> f32 {
        // OPLS-AA parameters (kJ/mol)
        match self {
            Self::Oxygen => 0.71128,
            Self::Hydrogen => 0.19246,
            Self::Carbon => 0.27614,
            Self::Nitrogen => 0.71128,
        }
    }
    pub fn sigma(&self) -> f32 {
        match self {
            Self::Oxygen => 0.3166,
            Self::Hydrogen => 0.225,
            Self::Carbon => 0.350,
            Self::Nitrogen => 0.325,
        }
    }
}

#[derive(Component, Default, Debug, Clone, Copy)]
pub struct Force {
    pub total: Vec3,
    pub bond: Vec3,
    pub angle: Vec3,
    pub non_bonded: Vec3,
}

impl Force {
    pub fn total_magnitude(&self) -> f32 {
        self.total.length()
    }
}

#[derive(Component, Deref, DerefMut, Debug)]
pub struct Velocity(pub Vec3);

#[derive(Component, Deref, DerefMut)]
pub struct Acceleration(pub Vec3);

// --- Components related to visualization ---

#[derive(Component)]
pub struct BondVisualization {
    pub atom1: Entity,
    pub atom2: Entity,
}
