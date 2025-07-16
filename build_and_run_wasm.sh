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

echo "--- Creating index.html ---"
# This creates the necessary index.html file in the static directory
cat > static/index.html <<- EOM
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8" />
    <title>Molecular Dynamics Simulation</title>
    <style>
      body { margin: 0; background-color: black; }
      canvas#bevy { width: 100vw; height: 100vh; }
    </style>
  </head>
  <body>
    <canvas id="bevy"></canvas>
    <script>
    const SAVE_TIMESTAMP_KEY = "simulation_save_timestamp";
    function js_save_to_localstorage(key, value) {
    console.log(`JS: Saving state to localStorage with key "${key}"`);
    window.localStorage.setItem(key, value);
    window.localStorage.setItem(SAVE_TIMESTAMP_KEY, new Date().toISOString());
    }
    function js_load_from_localstorage(key) {
    console.log(`JS: Loading state from localStorage with key "${key}"`);
    return window.localStorage.getItem(key);
    }
    function js_get_save_timestamp() {
    return window.localStorage.getItem(SAVE_TIMESTAMP_KEY);
    }
    </script>
    <script>
      const canvas = document.getElementById('bevy');
      canvas.addEventListener('contextmenu', (event) => {
        event.preventDefault();
      });
    </script>
    <script type="module">
      import init from './molecular_dynamics.js';
      init();
    </script>
  </body>
</html>
EOM

echo "--- Starting local web server on http://localhost:8080 ---"
# Serve the contents of the 'static' directory
python3 -m http.server 8080 --directory static
