use super::resources::*;
use crate::spawning::SpawnMoleculeFromSMILESEvent;
use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCamera;
pub struct SetupPlugin;

impl Plugin for SetupPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(
            Startup,
            (
                setup_scene,
                setup_shared_assets,
                trigger_initial_molecule_spawn,
            )
                .chain(),
        );
    }
}

fn setup_scene(mut commands: Commands) {
    info!("Setting up scene: Camera, Light, UI");
    commands.spawn((
        PanOrbitCamera::default(),
        Transform::from_xyz(-1.0, 2.0, 2.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        DirectionalLight::default(),
        Transform::from_xyz(10.0, 20.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
}

fn setup_shared_assets(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let bond_mesh_handle = meshes.add(Cylinder::new(0.015, 1.0));
    let bond_material_handle = materials.add(Color::srgb(0.8, 0.8, 0.2));
    let undefined_bond_material_handle = materials.add(Color::srgb(0.8, 0.2, 0.2)); // A bright red

    commands.insert_resource(SharedAssetHandles {
        bond_mesh: bond_mesh_handle,
        bond_material: bond_material_handle,
        undefined_bond_material: undefined_bond_material_handle,
    });
}

fn trigger_initial_molecule_spawn(
    mut spawn_writer: EventWriter<SpawnMoleculeFromSMILESEvent>,
    molecule_selection: Res<MoleculeSelection>,
) {
    let smiles = match molecule_selection.0.as_str() {
        "ethanol" => "CCO",
        "water" => "O",
        "alanine" => "C[C@H](N)C(=O)O",
        "glycine" => "NCC(=O)O",
        "acetaldehyde" => "CC=O",
        "arginine" => "N[C@@H](CCCNC(N)=N)C(O)=O",
        "phenylalanine" => "N[C@@H](CC1=CC=CC=C1)C(=O)O",
        "methylamine" => "CN",
        _ => "CCO", // Default to ethanol
    };
    info!(
        "Sending initial spawn event for '{}' with SMILES: {}",
        molecule_selection.0, smiles
    );
    spawn_writer.write(SpawnMoleculeFromSMILESEvent(smiles.to_string()));
}
