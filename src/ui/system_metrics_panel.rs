// src/ui/system_metrics_panel.rs

use crate::resources::EnergyHistory;
use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};
use egui_plot::{Legend, Line, Plot, PlotPoints};

pub struct SystemMetricsPanelPlugin;

impl Plugin for SystemMetricsPanelPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(EguiPrimaryContextPass, system_metrics_panel_egui_system);
    }
}

fn system_metrics_panel_egui_system(mut contexts: EguiContexts, history: Res<EnergyHistory>) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };

    egui::Window::new("System Metrics")
        // THE FIX: Change the anchor offset.
        // The "Selection Info" panel has a width of 250. We add some padding.
        // This places the metrics panel to the right of the info panel.
        .anchor(egui::Align2::LEFT_CENTER, egui::vec2(10.0, 40.0))
        .default_width(350.0)
        .default_height(200.0)
        .show(ctx, |ui| {
            let potential_line = Line::new(
                "Potential",
                PlotPoints::from_iter(history.potential.iter().map(|(x, y)| [*x, *y])),
            );
            let kinetic_line = Line::new(
                "Kinetic",
                PlotPoints::from_iter(history.kinetic.iter().map(|(x, y)| [*x, *y])),
            );
            let total_line = Line::new(
                "Total",
                PlotPoints::from_iter(history.total.iter().map(|(x, y)| [*x, *y])),
            );

            Plot::new("energy_plot")
                .legend(Legend::default())
                .x_axis_label("Simulation Step")
                .y_axis_label("Energy")
                .show(ui, |plot_ui| {
                    plot_ui.line(potential_line);
                    plot_ui.line(kinetic_line);
                    plot_ui.line(total_line);
                });
        });
}
