# CALPHAD server

Load any `.tdb` thermodynamic database once, then ask it questions over HTTP/JSON
from Python, C++, or any other language. pycalphad does the thermodynamics.

| File | What it is |
|---|---|
| `calphad_core.py` | All the thermodynamics. Importable on its own (no server needed). |
| `calphad_kks.py` | KKS phase-field building blocks for a binary, two-phase model (see [KKS](#kks-phase-field-building-blocks)). |
| `calphad_server.py` | HTTP server around the core (FastAPI). |
| `calphad_client.py` | Python client, standard library only. Copy into any project. |
| `cpp/` | C++ client (`calphad_client.hpp`), KKS header (`kks_parabolic.hpp`), two examples, `CMakeLists.txt`. |
| `examples/kks_1d_directional.py` | 1-D Ni–10Al directional solidification in a moving temperature gradient. |
| `selftest.py` | Consistency checks; run it on any new database. |
| `learning/` | One runnable script per part of the study guide (see below). |
| `cuni_cost507.tdb` | Cu–Ni (LIQUID, FCC) from COST 507. |
| `alni_dupin_2001.tdb` | Al–Ni (Dupin et al. 2001), from pycalphad's test databases. |
| `alcrni.tdb` | Al–Cr–Ni from pycalphad's test databases. Useful for exercising ordered phases (L1₂, B2); not a database for alloy design. |

## Learning materials

* **CALPHAD Study Guide: From Gibbs Energy to the KKS Server** (Claude Doc, https://claude.ai/artifact/RqGi3rQfcCsymzPJJ1DiAa): nine parts from Gibbs energy to reading this code, with exercises and answers. The scripts in `learning/` reproduce every snippet; run them from any folder, e.g. `python learning/part5_pycalphad.py`.
* **CALPHAD Server & KKS Framework Guide** (Claude Doc, https://claude.ai/artifact/DYaPx84Vvd9nVVGhyTVKLF): the full reference for every API, option and limit.

## Install and run

```bash
pip install pycalphad fastapi uvicorn matplotlib   # matplotlib: demo figures only
python calphad_server.py --tdb cuni_cost507.tdb --tdb alni_dupin_2001.tdb   # port 8000
python selftest.py                                                  # optional checks
```

Open <http://127.0.0.1:8000/docs> for interactive docs where you can try every call.
The server listens only on your machine by default; `--host 0.0.0.0` exposes it to
your network (there is no authentication, so only do that on a trusted network).

## Calling it

**Python**
```python
from calphad_client import CalphadClient
c = CalphadClient("http://127.0.0.1:8000", db="cuni_cost507")
r = c.equilibrium({"NI": 0.33}, T=1500, balance="CU")
for ph in r["results"]["phases"]:
    print(ph["phase"], ph["fraction"], ph["composition"])
```

**C++** (downloads two header-only libraries on first build)
```bash
cd cpp && cmake -S . -B build && cmake --build build && ./build/calphad_example
```
```cpp
calphad::CalphadClient c("127.0.0.1", 8000, "cuni_cost507");
auto r = c.equilibrium({{"NI", 0.33}}, 1500.0, "CU");
double fL = r["results"]["phases"][0]["fraction"];
```

**Anything else**: send JSON with any HTTP library, e.g.
```bash
curl -X POST http://127.0.0.1:8000/equilibrium -H 'Content-Type: application/json' \
     -d '{"db":"cuni_cost507","composition":{"NI":0.33},"balance":"CU","T":1500}'
```
`/openapi.json` describes every endpoint; tools such as openapi-generator turn it
into typed clients for Java, C#, Julia, Rust, etc.

## Endpoints

| Method | Path | Purpose |
|---|---|---|
| GET | `/health` | Server up? Which databases? |
| GET | `/databases` | Loaded databases with elements and phases |
| POST | `/databases` | Load another TDB: `{"path": "..."}` (server-side file) or `{"tdb_text": "...", "name": "..."}` |
| GET | `/databases/{db}` | Elements, atomic masses, phases, sublattice models, order/disorder pairs |
| POST | `/equilibrium` | **API 1** |
| POST | `/phase_energies` | **API 2** |
| POST | `/energy_surface` | **API 3** |
| POST | `/kks/parameters` | KKS inputs: tie lines, parabolic coefficients, G tables (see below) |

`db` may be left out of a request when only one database is loaded.

### Specifying the alloy (APIs 1 and 2)

```json
{"composition": {"AL": 0.2, "CR": 0.1}, "balance": "NI"}        // solutes + balance
{"composition": {"CU": 0.67, "NI": 0.33}}                       // all elements, sum 1
{"composition": {"NI": 30}, "balance": "CU", "units": "weight_percent"}
```
`units`: `mole_fraction` (default), `atomic_percent`, `weight_fraction`, `weight_percent`.
Only the elements you name are included in the calculation, so a large
database can be used for any subsystem. Optional: `"P"` (Pa, default 101325),
`"phases"` (restrict the candidate phases).

### API 1 — `POST /equilibrium`

Stable phases, their fractions and compositions at fixed alloy composition.
`T` can be a number or a list (one result per temperature).

```json
{"db": "alcrni", "composition": {"AL": 0.2, "CR": 0.1}, "balance": "NI", "T": 1600}
```
Response (`results` is a list when `T` was a list):
```json
{"alloy": {"AL": 0.2, "CR": 0.1, "NI": 0.7},
 "results": {
   "T": 1600.0, "GM_system": -118173.48,
   "chemical_potentials": {"AL": -189065.45, "CR": -93171.75, "NI": -101490.30},
   "phases": [
     {"phase": "L12_FCC", "composition_set": 1, "state": "disordered", "equivalent_to": "FCC_A1",
      "fraction": 0.8555, "composition": {"AL": 0.1857, "CR": 0.1073, "NI": 0.7070},
      "composition_wt": {"...": "..."}, "GM": -116858.1,
      "site_fractions": [{"AL": 0.186, "CR": 0.107, "NI": 0.707},
                         {"AL": 0.186, "CR": 0.107, "NI": 0.707}]},
     {"phase": "B2", "composition_set": 1, "state": "ordered", "equivalent_to": "B2",
      "fraction": 0.1445, "composition": {"AL": 0.2848, "CR": 0.0567, "NI": 0.6585},
      "GM": -125959.1, "...": "..."}],
   "stability_check": {"verified": true, "sampling_density": 2000,
                       "driving_forces_absent_phases": {"LIQUID": -26.09}}}}
```
(Numbers from the bundled `alcrni.tdb`; values abbreviated.)
* `fraction`: fraction of the total moles of atoms in that phase.
* `GM` per phase: its molar Gibbs energy at its own composition (= Σ xᵢμᵢ).
* The same phase can appear twice (`composition_set` 2): a miscibility gap, or
  ordered + disordered states of one order/disorder phase. `state` /
  `equivalent_to` say which, e.g. `L12_FCC` that is disordered is really γ (`FCC_A1`).

### API 2 — `POST /phase_energies`

Gibbs energy of **every** candidate phase at the alloy composition, stable or not.
Each phase is minimised on its own, as one homogeneous phase (its sublattice
occupancies and ordering relax; its composition is held at the alloy's).

```json
{"db": "cuni_cost507", "composition": {"NI": 0.33}, "balance": "CU", "T": 1500}
```
Per phase:
* `GM` — J/mol of atoms.
* `chemical_potentials` — this phase's tangent plane at this composition.
  `μ_i − μ_dependent` are the diffusion potentials a KKS solver equates between phases.
* `GM_above_equilibrium` — how far above the equilibrium mixture (≥ 0).
* `driving_force` — against the equilibrium chemical potentials; 0 for stable
  phases, negative for phases that would not form.
* `splits_if_allowed` — present when the phase alone would rather separate into
  two compositions (miscibility gap, or ordered + disordered).

### API 3 — `POST /energy_surface`

G(x) of every phase at temperature `T`: curves for a binary, surfaces for a ternary.

```json
{"db": "alcrni", "T": 1273, "step": 0.02}
{"db": "alcrni", "T": 1273, "elements": ["AL", "NI"]}         // a binary subsystem
{"db": "alcrni", "T": 1273, "mode": "points", "pdens": 2000}  // raw sampled configurations
```
Response per phase: `X` (rows of mole fractions, columns in the order of
`elements`), `GM` (J/mol), `approximate`.
* `mode: "grid"` (default): values on a regular composition grid, spacing `step`.
  Exact for substitutional solutions. For phases with sublattices it is the lowest
  sampled G in each grid cell (`"approximate": true`; raise `pdens` to refine).
* `mode: "points"`: every sampled internal configuration (more points, includes
  higher-energy orderings).

## Conventions

* Energies: J per mole of atoms, relative to the database's reference state
  (normally SER). Only differences between phases/compositions are physical.
* Compositions returned are mole fractions; `composition_wt` gives mass fractions.
* T in K, P in Pa.

## Reliability: the stability check

pycalphad starts its minimiser from a sampled set of configurations. If the
sample misses the low-energy region of a phase, it can return a local minimum.
This happened in testing: Ni–20Al–10Cr at 1600 K came back as LIQUID + L1₂, while
L1₂ + B2 is lower in energy. So every equilibrium is checked: each absent phase's
driving force must be ≤ 1 J/mol. If one is positive, the calculation is repeated
with denser sampling (2000, then 10000 points). `stability_check` reports what
happened; `"verified": false` means even the densest sampling left a phase with a
positive driving force, so treat that result with care. Pass `"verify": false`
to skip the check when speed matters more.

## Using it from a phase-field / KKS code

An HTTP round trip costs about a millisecond, which is fine for setup but too slow
to call at every grid point every time step. Typical pattern:

1. **Before the run**, fetch what the solver needs: tie-line compositions vs T
   (API 1 with a list of T), G and diffusion potentials at chosen compositions
   (API 2), or G(x) tables to fit or interpolate (API 3).
2. **In the time loop**, use those fits/tables (parabolic or spline free energies).

For direct in-process use from Python, skip the server:
```python
from calphad_core import CalphadSystem
s = CalphadSystem.from_file("cuni_cost507.tdb")
s.equilibrium({"NI": 0.33}, 1500, balance="CU")
```

## KKS phase-field building blocks

`calphad_kks.py` turns a TDB into what a KKS solver needs, for a binary alloy
and two phases (e.g. γ = `FCC_A1` and `LIQUID`). The server's three APIs are
unaffected; this is a separate module plus one extra endpoint.

```python
from calphad_kks import KKSModel
kks = KKSModel.from_tdb("alni_dupin_2001.tdb", solute="AL",
                        phase_a="FCC_A1", phase_b="LIQUID",   # phi = 1 is phase_a
                        mode="exact", Vm=7.0e-6, c_alloy=0.10)

# in the time loop; c, phi, T are NumPy arrays (T can be a field)
ca, cb, mu = kks.split(c, phi=phi, T=T, ca0=ca, cb0=cb)   # warm start from last step
dF = kks.driving_force(ca, cb, mu, T)                    # J/m^3
# dphi/dt = M [ eps2 lap(phi) - W dg(phi) - dh(phi) dF ]
# dc/dt   = div[ D(phi) ( h grad(ca) + (1-h) grad(cb) ) ]
eps2, W = KKSModel.interface_parameters(sigma=0.3, width=25e-9)
```

| Method | Returns |
|---|---|
| `split(c, phi, T)` | c_a, c_b, diffusion potential μ (J/mol): mixture rule + equal μ at every point |
| `driving_force(ca, cb, mu, T)` | f_a − f_b − μ(c_a − c_b), J/m³ (negative favours phase a) |
| `dmu_dc(ca, cb, phi, T)` | dμ/dc, for the μ-form of the diffusion equation |
| `tie_lines(T_list)` | c_a, c_b, k, μ, G and curvatures A at each T, checked against pycalphad (`consistent`) |
| `h, dh, g, dg` | interpolation and double-well functions |
| `interface_parameters(σ, width)` | ε², W from interface energy and width (KKS 1999 relations) |
| `export_parabolic_csv`, `export_json`, `parameters` | the same data for other codes |

**Free-energy modes**

| `mode` | What it is | Use when |
|---|---|---|
| `exact` | the CALPHAD expression compiled to machine code, exact derivatives | substitutional phases (liquid, γ, …); most accurate |
| `table` | G(c, T) computed from CALPHAD and splined | any phase, including ordered ones (γ′ etc.) |
| `parabolic` | G_eq + μ_eq(c − c_eq) + A/2(c − c_eq)² about the tie line at each T | speed; exact at the tie line, approximate away from it |

`table` and `parabolic` need a temperature list, e.g. `T=(1690, 1725, 1)`.
For order/disorder phases pass `branches={"FCC_L12": "ordered"}` to tabulate the
ordered branch only (γ′); otherwise the lower of the two states is used and
`kks.warnings` tells you G has a kink.

**Speed** (512 × 512 grid with a temperature field, this machine): `split` takes
0.22 s cold / 0.06 s warm-started (`exact`), 1.1 s (`table`), 0.03 s (`parabolic`).

**From C++**: fetch coefficients once with `POST /kks/parameters` (or
`KKSModel.export_parabolic_csv`), then use `cpp/kks_parabolic.hpp`, which does
the split, driving force and dμ/dc in closed form with T interpolation.
`cpp/kks_example.cpp` shows both steps; its numbers match the Python module.

**Demo**: `python examples/kks_1d_directional.py` runs planar directional
solidification of Ni–10Al at G = 10⁷ K/m, V = 0.01 m/s. Solute is conserved to
~10⁻¹³; after the classic initial transient the interface settles 0.5 K below
the solidus of the alloy (kinetic undercooling from the finite interface
mobility), the solid approaches c0 and the interface liquid approaches c0/k.

Validated against pycalphad: KKS tie lines match equilibrium to 10⁻⁵ or better
(Cu–Ni, γ/L and γ′/L in Ni–Al); the split satisfies mass balance to 10⁻¹⁶ and
equal potentials to 10⁻⁸ J/mol; tabulated energies match exact ones to 10⁻³ J/mol.

**Not covered**: more than two phases, more than two elements, and
anti-trapping currents / thin-interface calibration of M (those belong in your
solver; see the notes in `calphad_kks.py`).

## Limitations

* The stability check uses the larger of pycalphad's driving force and a sampled
  estimate for each absent phase. It does not catch every missed miscibility gap
  inside a phase that is already stable.
* API 2 refines each phase from several ordered and disordered starting points,
  because pycalphad alone can return the metastable disordered state of an
  ordered phase (Ni–25Al γ′ at 1600 K came back 355 J/mol too high).
* Ionic-liquid and other exotic models work as far as pycalphad supports them;
  only the substitutional and sublattice (CEF) models were tested here.
* No authentication: keep the default `127.0.0.1` binding unless you trust the network.
