// src/ui/amino_acid_spawner.rs

use crate::resources::{AminoAcidFile, AminoAcidInfo};
use crate::spawning::SpawnMoleculeFromSMILESEvent;
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

pub struct AminoAcidSpawnerPlugin;

#[derive(Resource, Default)]
struct SpawnerState {
    is_open: bool,
}

impl Plugin for AminoAcidSpawnerPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<SpawnerState>()
            .add_systems(Startup, setup_amino_acid_info)
            .add_systems(
                EguiPrimaryContextPass,
                (
                    amino_acid_spawner_egui_system,
                    amino_acid_spawner_toggle_button,
                ),
            );
    }
}

// On startup, load the amino acid definitions from the JSON file.
fn setup_amino_acid_info(mut commands: Commands) {
    #[cfg(target_arch = "wasm32")]
    let file_contents = include_str!("../../assets/amino_acids.json");

    #[cfg(not(target_arch = "wasm32"))]
    let file_contents = std::fs::read_to_string("assets/amino_acids.json")
        .expect("Failed to read amino_acids.json");

    // 1. Deserialize into the struct that matches the file's structure.
    let aa_file: AminoAcidFile =
        serde_json::from_str(&file_contents).expect("Failed to parse amino_acids.json");

    // 2. Insert the inner vector into the AminoAcidInfo resource.
    commands.insert_resource(AminoAcidInfo(aa_file.amino_acids));
}
fn amino_acid_spawner_egui_system(
    mut contexts: EguiContexts,
    mut state: ResMut<SpawnerState>,
    amino_acid_info: Res<AminoAcidInfo>,
    mut spawn_writer: EventWriter<SpawnMoleculeFromSMILESEvent>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    egui::Window::new("Amino Acid Spawner")
        .open(&mut state.is_open)
        .default_pos(ctx.screen_rect().right_bottom() + egui::vec2(10.0, 10.0))
        .show(ctx, |ui| {
            ui.heading("Standard Amino Acids");
            ui.label("Click to spawn a molecule.");
            ui.add_space(5.0);

            // Create a scrolling area to fit all the buttons.
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.set_width(200.0);
                // Create buttons for each amino acid in a grid.
                egui::Grid::new("amino_acid_grid").show(ui, |ui| {
                    for (i, acid) in amino_acid_info.0.iter().enumerate() {
                        if ui
                            .button(format!("{} ({})", acid.name, acid.three_letter_code))
                            .clicked()
                        {
                            // When clicked, send the existing event with the correct SMILES.
                            spawn_writer.write(SpawnMoleculeFromSMILESEvent(
                                acid.smiles_zwitterionic.clone(),
                                None, // `None` specifies to spawn just one molecule.
                            ));
                        }
                        if (i + 1) % 2 == 0 {
                            // 2 buttons per row
                            ui.end_row();
                        }
                    }
                });
            });
        });
}

fn amino_acid_spawner_toggle_button(mut contexts: EguiContexts, mut state: ResMut<SpawnerState>) {
    let Ok(ctx) = contexts.ctx_mut() else { return };
    egui::Area::new(egui::Id::new("spawner_toggle_button_area"))
        .anchor(egui::Align2::RIGHT_CENTER, egui::vec2(10.0, 10.0))
        .show(ctx, |ui| {
            if ui.button("Amino Acids").clicked() {
                state.is_open = !state.is_open;
            }
        });
}
