use bevy::prelude::*;

#[derive(Component, Debug, Clone)]
pub struct Atom {
    pub type_name: String,
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

#[derive(Component)]
pub struct BondVisualization {
    pub atom1: Entity,
    pub atom2: Entity,
    pub strand_index: u8,
    pub total_strands: u8,
}
