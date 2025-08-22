from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware

# What the client sends
class MoleculeData(BaseModel):
    smiles: str
    mol: str

app = FastAPI(title="Molecule API", version="0.1.0")

# Allow the client origins (adjust as needed)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://127.0.0.1:8000", "http://localhost:8000"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def health():
    return {"ok": True}

@app.post("/api/analyze")
async def analyze_molecule(data: MoleculeData):
    # Minimal echo back (replace with your real logic later)
    print("Got SMILES:", data.smiles)  # 👈 just prints it in the backend terminal
    return {
        "message": "Molecule received successfully!",
        "smiles": data.smiles,
        "mol_length": len(data.mol)
    }