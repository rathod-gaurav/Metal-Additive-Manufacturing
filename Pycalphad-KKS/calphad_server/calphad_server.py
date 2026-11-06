"""
calphad_server.py
=================
HTTP/JSON server around calphad_core, so any language can query a CALPHAD
database.

Start
-----
  python calphad_server.py --tdb cuni_cost507.tdb --tdb alcrni.tdb
  python calphad_server.py --tdb mydb.tdb --host 0.0.0.0 --port 8000

Then open http://127.0.0.1:8000/docs for interactive documentation; the
machine-readable spec at /openapi.json can generate clients for C++, Java,
Julia, MATLAB, etc. (e.g. with openapi-generator).

Endpoints
---------
  GET  /health
  GET  /databases                  loaded databases
  POST /databases                  load another TDB (by path or text)
  GET  /databases/{db}             elements, phases, sublattice models
  POST /equilibrium                API 1: phase fractions + compositions
  POST /phase_energies             API 2: G of every phase at a composition
  POST /energy_surface             API 3: G(x) surfaces for all phases
  POST /kks/parameters             KKS inputs for a binary, two-phase model
                                   (tie lines, parabolic coefficients, G tables)

Every POST body may name the database with "db"; if only one database is
loaded it is used by default.
"""

import argparse
import os
import threading
from typing import Dict, List, Literal, Optional, Union

import numpy as np
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from calphad_core import P_DEFAULT, CalphadError, CalphadSystem
from calphad_kks import KKSModel

app = FastAPI(
    title="CALPHAD server",
    version="1.0",
    description="Equilibrium, phase Gibbs energies and Gibbs energy surfaces "
                "from any TDB database (pycalphad underneath). Energies in J/mol "
                "of atoms, compositions in mole fractions, T in K, P in Pa.",
)

_DBS: Dict[str, CalphadSystem] = {}
_DBS_LOCK = threading.Lock()


def _get_db(name: Optional[str]) -> CalphadSystem:
    if name is None:
        if len(_DBS) == 1:
            return next(iter(_DBS.values()))
        raise HTTPException(400, f"Several databases loaded {sorted(_DBS)}; say which with 'db'.")
    if name not in _DBS:
        raise HTTPException(404, f"Database '{name}' not loaded. Loaded: {sorted(_DBS)}")
    return _DBS[name]


def _run(fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except CalphadError as exc:
        raise HTTPException(400, str(exc))
    except HTTPException:
        raise
    except Exception as exc:  # surface pycalphad errors to the caller
        raise HTTPException(500, f"{type(exc).__name__}: {exc}")


# ---------------------------------------------------------------------------
# Request models
# ---------------------------------------------------------------------------
Units = Literal["mole_fraction", "atomic_percent", "weight_fraction", "weight_percent"]


class AlloyRequest(BaseModel):
    db: Optional[str] = Field(None, description="Database name (see GET /databases)")
    composition: Dict[str, float] = Field(
        ..., examples=[{"NI": 0.3, "CU": 0.7}],
        description="Element amounts. Either all elements (summing to 1 or 100), "
                    "or the solutes only together with 'balance'.")
    units: Units = "mole_fraction"
    balance: Optional[str] = Field(None, description="Element that makes up the rest, e.g. 'NI'")
    P: float = P_DEFAULT
    phases: Optional[List[str]] = Field(None, description="Restrict to these phases (default: all)")


class EquilibriumRequest(AlloyRequest):
    T: Union[float, List[float]] = Field(..., examples=[1500.0],
                                         description="K. A list returns one result per T.")
    verify: bool = Field(True, description="Check no absent phase has a positive driving "
                                           "force; re-solve with denser sampling if one does.")


class PhaseEnergiesRequest(AlloyRequest):
    T: float = Field(..., examples=[1500.0])


class SurfaceRequest(BaseModel):
    db: Optional[str] = None
    T: float = Field(..., examples=[1500.0])
    P: float = P_DEFAULT
    elements: Optional[List[str]] = Field(None, description="Subsystem, e.g. ['CU','NI'] (default: all)")
    phases: Optional[List[str]] = None
    mode: Literal["grid", "points"] = "grid"
    step: float = Field(0.02, gt=0, le=0.5, description="Grid spacing in mole fraction")
    pdens: int = Field(500, gt=0, description="Sampling density for phases with sublattices")
    max_points: int = Field(200000, gt=0)


class LoadRequest(BaseModel):
    name: Optional[str] = Field(None, description="Name to register it under")
    path: Optional[str] = Field(None, description="TDB path on the server's machine")
    tdb_text: Optional[str] = Field(None, description="Full TDB file contents")


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health():
    return {"status": "ok", "databases": sorted(_DBS)}


@app.get("/databases")
def list_databases():
    return {name: {"elements": s.elements, "phases": s.phases} for name, s in _DBS.items()}


@app.post("/databases")
def load_database(req: LoadRequest):
    if bool(req.path) == bool(req.tdb_text):
        raise HTTPException(400, "Give exactly one of 'path' or 'tdb_text'.")
    if req.path:
        system = _run(CalphadSystem.from_file, req.path, req.name)
    else:
        system = _run(CalphadSystem.from_text, req.tdb_text, req.name or "db")
    with _DBS_LOCK:
        _DBS[system.name] = system
    return system.info()


@app.get("/databases/{name}")
def database_info(name: str):
    return _get_db(name).info()


@app.post("/equilibrium")
def api_equilibrium(req: EquilibriumRequest):
    """API 1: stable phases, phase fractions and phase compositions."""
    s = _get_db(req.db)
    return _run(s.equilibrium, req.composition, req.T, req.P, req.phases, req.units,
                req.balance, req.verify)


@app.post("/phase_energies")
def api_phase_energies(req: PhaseEnergiesRequest):
    """API 2: Gibbs energy of every phase at the alloy composition."""
    s = _get_db(req.db)
    return _run(s.phase_energies, req.composition, req.T, req.P, req.phases, req.units, req.balance)


@app.post("/energy_surface")
def api_energy_surface(req: SurfaceRequest):
    """API 3: G(x) for every phase over the composition space at T."""
    s = _get_db(req.db)
    return _run(s.energy_surface, req.T, req.P, req.elements, req.phases,
                req.step, req.mode, req.pdens, req.max_points)


class KKSRequest(BaseModel):
    db: Optional[str] = None
    solute: str = Field(..., examples=["AL"], description="c = mole fraction of this element")
    solvent: Optional[str] = Field(None, description="Other element (needed if the database has more than two)")
    phase_a: str = Field(..., examples=["FCC_A1"], description="Phase at phi = 1")
    phase_b: str = Field(..., examples=["LIQUID"], description="Phase at phi = 0")
    T: List[float] = Field(..., examples=[[1700.0, 1705.0, 1710.0]],
                           description="Temperatures for tie lines / tables (K)")
    mode: Literal["exact", "table", "parabolic"] = "parabolic"
    c_alloy: Optional[float] = Field(None, description="Nominal alloy composition; picks the "
                                                       "right two-phase field")
    Vm: Optional[float] = Field(None, description="Molar volume, m^3/mol")
    branches: Optional[Dict[str, Literal["lowest", "ordered", "disordered"]]] = None
    include_tables: bool = Field(False, description="Also return G, dG/dc, d2G/dc2 on a c grid")
    n_c: int = Field(199, gt=7, le=2001, description="Points in the returned c grid")


_KKS_CACHE: Dict[tuple, KKSModel] = {}


@app.post("/kks/parameters")
def api_kks_parameters(req: KKSRequest):
    """KKS inputs for a binary two-phase model: tie lines vs T (c_a, c_b, mu, G,
    curvatures A = d2G/dc2, partition coefficient), and optionally tabulated
    G(c), dG/dc, d2G/dc2 for both phases. Fetch once at solver start-up."""
    s = _get_db(req.db)
    temps = sorted(float(t) for t in req.T)
    key = (s.name, req.solute.upper(), (req.solvent or "").upper(), req.phase_a.upper(),
           req.phase_b.upper(), req.mode, req.c_alloy, req.Vm,
           tuple(sorted((req.branches or {}).items())),
           None if req.mode == "exact" else tuple(temps))
    if key not in _KKS_CACHE:
        _KKS_CACHE[key] = _run(KKSModel, s, req.solute, req.phase_a, req.phase_b, mode=req.mode,
                               Vm=req.Vm, c_alloy=req.c_alloy,
                               T=None if req.mode == "exact" else temps,
                               solvent=req.solvent, branches=req.branches)
    kks = _KKS_CACHE[key]
    c_grid = np.linspace(1e-3, 1 - 1e-3, req.n_c)
    out = _run(kks.parameters, temps, c_grid, req.include_tables)
    out["warnings"] = kks.warnings
    return out


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="CALPHAD HTTP server")
    ap.add_argument("--tdb", action="append", default=[],
                    help="TDB file to load (repeatable). Registered under its file name.")
    ap.add_argument("--host", default="127.0.0.1",
                    help="127.0.0.1 = this machine only; 0.0.0.0 = whole network")
    ap.add_argument("--port", type=int, default=8000)
    args = ap.parse_args()
    for path in args.tdb:
        s = CalphadSystem.from_file(path)
        _DBS[s.name] = s
        print(f"Loaded {path} as '{s.name}': {s.elements}, {len(s.phases)} phases")
    import uvicorn
    uvicorn.run(app, host=args.host, port=args.port, log_level="info")


if __name__ == "__main__":
    main()
