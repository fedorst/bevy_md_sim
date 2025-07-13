// src/ui/peptide_builder.rs

use crate::resources::{AminoAcidFile, AminoAcidInfo};
use crate::spawning::SpawnMoleculeFromSMILESEvent;
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};

pub struct AminoAcidSpawnerPlugin;

impl Plugin for AminoAcidSpawnerPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_amino_acid_info)
            .add_systems(EguiPrimaryContextPass, amino_acid_spawner_egui_system);
    }
}

// On startup, load the amino acid definitions from the JSON file.
fn setup_amino_acid_info(mut commands: Commands) {
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
    amino_acid_info: Res<AminoAcidInfo>,
    // We now send the standard spawn event directly.
    mut spawn_writer: EventWriter<SpawnMoleculeFromSMILESEvent>,
) {
    let Ok(ctx) = contexts.ctx_mut() else { return };

    egui::Window::new("Amino Acid Spawner")
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
