use serde::Deserialize;

#[derive(Deserialize, Debug)]
pub struct AtomSpec {
    pub id: String,
    pub type_name: String,
    pub pos: [f32; 3],
    pub element: String,
}

#[derive(Deserialize, Debug)]
pub struct BondSpec {
    pub atoms: [String; 2],
    pub order: String,
}

#[derive(Deserialize, Debug)]
pub struct MoleculeConfig {
    pub name: String,
    pub atoms: Vec<AtomSpec>,
    pub bonds: Vec<BondSpec>,
}
