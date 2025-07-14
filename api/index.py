from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
# Assuming auto_typer.py is also in the api/ directory
from auto_typer import build_molecule_from_smiles

app = FastAPI()

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], # Start with a wildcard for testing
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class SmilesRequest(BaseModel):
    smiles: str

@app.post("/api/generate_molecule") # The full URL will be your-site.vercel.app/api/generate_molecule
def generate_molecule(request: SmilesRequest):
    # The endpoint now directly calls your existing, well-tested function
    return build_molecule_from_smiles(request.smiles, "Generated Molecule")
