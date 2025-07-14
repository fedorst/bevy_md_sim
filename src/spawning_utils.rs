// src/spawning_utils.rs

use bevy::prelude::*;

/// A pure function that returns the standard radius and color for a given element.
/// This is the single source of truth for atom visuals.
pub fn get_atom_visuals(element: &str) -> (f32, Color) {
    let radius = match element {
        "C" => 0.06,
        "O" => 0.05,
        "H" => 0.03,
        "N" => 0.055,
        "S" => 0.1,
        _ => 0.045,
    };
    let color = match element {
        "C" => Color::srgb(0.2, 0.2, 0.2),
        "O" => Color::srgb(1.0, 0.1, 0.1),
        "H" => Color::srgb(0.9, 0.9, 0.9),
        "N" => Color::srgb(0.1, 0.1, 1.0),
        "S" => Color::srgb(1.0, 1.0, 0.0),
        _ => Color::srgb(1.0, 0.2, 0.8),
    };
    (radius, color)
}
