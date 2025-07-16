# Bevy Molecular Dynamics
Real-time molecular dynamics simulation in Rust using the Bevy engine. Simulate molecules from SMILES strings with custom physics and interactive manipulation.

**Demo:** [md.fedor.ee](https://md.fedor.ee)
![Simulation Demo](https://github.com/fedorst/bevy_md_sim/blob/main/assets/demo_screencast.gif)

## Technical Overview
- **Physics**: Bond, angle, Lennard-Jones, and Coulomb force calculations
- **Molecule generation**: SMILES → 3D via Python/RDKit backend
- **Interaction**: Atom selection, force application, simulation pausing
- **UI**: egui panels for amino acid spawning and force visualization

## Architecture
- **Native core**: Rust/Bevy desktop application
- **Web deployment**: WASM frontend + FastAPI backend
- **Infrastructure**: Docker containers, Nginx reverse proxy, DigitalOcean VPS
- **CI/CD**: GitHub Actions for automated builds and deployments

## Local Execution
```bash
pip install -r api/requirements.txt
cargo run
```
