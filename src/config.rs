use crate::core::AtomType;
use serde::Deserialize;

#[derive(Deserialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[serde(rename_all = "PascalCase")]
pub enum ConfigAtomType {
    Oxygen,
    Hydrogen,
    Carbon,
    Nitrogen,
}

impl From<ConfigAtomType> for AtomType {
    fn from(config_type: ConfigAtomType) -> Self {
        match config_type {
            ConfigAtomType::Oxygen => AtomType::Oxygen,
            ConfigAtomType::Hydrogen => AtomType::Hydrogen,
            ConfigAtomType::Carbon => AtomType::Carbon,
            ConfigAtomType::Nitrogen => AtomType::Nitrogen,
        }
    }
}

#[derive(Deserialize, Debug)]
pub struct AtomSpec {
    pub id: String,
    #[serde(rename = "type")]
    pub atom_type: ConfigAtomType,
    pub pos: [f32; 3],
}

#[derive(Deserialize, Debug)]
pub struct MoleculeConfig {
    pub name: String,
    pub atoms: Vec<AtomSpec>,
    pub bonds: Vec<[String; 2]>,
}
