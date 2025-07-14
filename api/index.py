from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from src_py.auto_typer import build_molecule_from_smiles # Import from your library
from rdkit import rdBase

# Vercel will look for this 'app' variable
app = FastAPI()

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class SmilesRequest(BaseModel):
    smiles: str

# The path here is what the app sees *after* the Vercel rewrite
@app.post("/generate_molecule")
def generate_molecule(request: SmilesRequest):
    """
    Generates a typed molecule JSON from a SMILES string.
    """
    try:
        molecule_data = build_molecule_from_smiles(request.smiles, "Generated Molecule")
        if molecule_data is None or "error" in molecule_data:
            error_message = molecule_data.get("error", "Unknown error during molecule generation.")
            raise HTTPException(status_code=400, detail=error_message)
        return molecule_data
    except Exception as e:
        # Catch any other unexpected errors during processing
        raise HTTPException(status_code=500, detail=f"An internal error occurred: {str(e)}")

# Optional: Add a root endpoint for health checks
@app.get("/")
def read_root():
    return {"status": "ok", "message": "Molecule Generation API is running."}

@app.get("/test_rdkit")
def test_rdkit_version():
    """A simple endpoint to verify that RDKit can be imported and used."""
    try:
        version = rdBase.rdkitVersion
        return {"status": "ok", "rdkit_version": version}
    except Exception as e:
        return {"status": "error", "detail": str(e)}
