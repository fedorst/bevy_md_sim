// src/ui/force_inspector.rs

use crate::components::{Atom, Force};
use crate::interaction::SelectionState;
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};
use egui_plot::{Bar, BarChart, Legend, Plot};

pub struct ForceInspectorPlugin;

impl Plugin for ForceInspectorPlugin {
    fn build(&self, app: &mut App) {
        // This system will draw our new plot.
        app.add_systems(EguiPrimaryContextPass, force_inspector_egui_system);
    }
}

fn force_inspector_egui_system(
    mut contexts: EguiContexts,
    selection: Res<SelectionState>,
    atom_query: Query<&Force, With<Atom>>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };

    if selection.selected.len() != 1 {
        return;
    }

    let Ok(force) = atom_query.get(selection.selected[0]) else {
        return;
    };

    egui::Window::new("Force Inspector")
        .anchor(egui::Align2::RIGHT_TOP, egui::vec2(-10.0, 10.0))
        .show(ctx, |ui| {
            let bond_force = force.bond.length();
            let angle_force = force.angle.length();
            let non_bonded_force = force.non_bonded.length();

            // Define the colors clearly upfront
            let color_bond = egui::Color32::from_rgb(100, 200, 100);
            let color_angle = egui::Color32::from_rgb(100, 100, 200);
            let color_non_bonded = egui::Color32::from_rgb(200, 100, 100);

            // --- MODIFICATION: Set the color on the BarChart as well ---

            // Chart 1: Bond
            let chart_bond = BarChart::new(
                "Bond",
                vec![Bar::new(0.0, bond_force as f64).fill(color_bond)],
            )
            .width(0.7)
            .color(color_bond); // Set the color for the legend entry

            // Chart 2: Angle
            let chart_angle = BarChart::new(
                "Angle",
                vec![
                    Bar::new(0.0, angle_force as f64)
                        .base_offset(bond_force as f64)
                        .fill(color_angle),
                ],
            )
            .width(0.7)
            .color(color_angle); // Set the color for the legend entry

            // Chart 3: Non-Bonded
            let chart_non_bonded = BarChart::new(
                "Non-Bonded",
                vec![
                    Bar::new(0.0, non_bonded_force as f64)
                        .base_offset((bond_force + angle_force) as f64)
                        .fill(color_non_bonded),
                ],
            )
            .width(0.7)
            .color(color_non_bonded); // Set the color for the legend entry

            Plot::new("force_breakdown_plot")
                .legend(Legend::default())
                .width(200.0)
                .show_x(false)
                .show_y(true)
                .y_axis_label("Force Magnitude")
                .show(ui, |plot_ui| {
                    plot_ui.bar_chart(chart_bond);
                    plot_ui.bar_chart(chart_angle);
                    plot_ui.bar_chart(chart_non_bonded);
                });
        });
}
