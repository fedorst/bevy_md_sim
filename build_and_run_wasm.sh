#!/bin/bash
set -e # Exit immediately if a command fails

echo "--- Building WASM binary for release ---"
cargo build --release --target wasm32-unknown-unknown

echo "--- Running wasm-bindgen ---"
wasm-bindgen --out-dir static --target web \
  target/wasm32-unknown-unknown/release/molecular_dynamics.wasm

echo "--- Copying assets to static directory ---"
# This is the key step to fix the cursor errors
cp -r assets static/

echo "--- Copying index.html ---"
cp index.html static/

echo "--- Starting local web server on http://localhost:8080 ---"
# Serve the contents of the 'static' directory
python3 -m http.server 8080 --directory static
