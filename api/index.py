from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
# Assuming auto_typer.py is also in the api/ directory
from src_py.auto_typer import build_molecule_from_smiles
import uvicorn

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

@app.post("/generate_molecule")
def generate_molecule(request: SmilesRequest):
    return build_molecule_from_smiles(request.smiles, "Generated Molecule")

@app.get("/api/health")
def health_check():
    return {"status": "healthy"}

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
