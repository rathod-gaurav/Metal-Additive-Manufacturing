# Part 7: the framework layer by layer (run from inside the calphad_server folder)
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
from calphad_core import CalphadSystem

# --- 7.3 Loading a database -------------------------------------------------
s = CalphadSystem.from_file("alcrni.tdb")
print("name:", s.name)
print("elements:", s.elements)
print("phases:", s.phases)
print("L12_FCC entry of info():", s.info()["phases"]["L12_FCC"])

# --- 7.4 From your input to pycalphad's conditions ---------------------------
x = s._parse_composition({"AL": 10, "CR": 8}, units="weight_percent", balance="NI")
print("\nmole fractions:", {e: round(val, 4) for e, val in x.items()})
comps = s._components(x)
print("components:", comps)
print("conditions:", s._conditions(x, 1500, 101325))
print("candidate phases:", s._phase_list(comps, None))

# --- 7.6 API 1 ------------------------------------------------------------------
import time
cuni = CalphadSystem.from_file("cuni_cost507.tdb")
t0 = time.time()
r = cuni.equilibrium({"NI": 0.33}, T=1500, balance="CU")
print(f"\nAPI 1 took {time.time()-t0:.1f} s")
res = r["results"]
print("GM_system:", round(res["GM_system"], 2), " mu:", {k: round(v, 1) for k, v in res["chemical_potentials"].items()})
for ph in res["phases"]:
    print(ph["phase"], "set", ph["composition_set"], "fraction", round(ph["fraction"], 4),
          "x_Ni", round(ph["composition"]["NI"], 4), "GM", round(ph["GM"], 2),
          "state", ph["state"])
print("stability_check:", res["stability_check"])

# --- 7.7 API 2 ------------------------------------------------------------------
alni = CalphadSystem.from_file("alni_dupin_2001.tdb")
t0 = time.time()
r2 = alni.phase_energies({"AL": 0.25}, T=1600, balance="NI")
print(f"\nAPI 2 took {time.time()-t0:.1f} s")
print("stable:", r2["stable_phases"], " GM_equilibrium:", round(r2["GM_equilibrium"], 2))
for p in r2["phases"]:
    if p["GM"] is None:
        print(f"{p['phase']:9s} GM None  ({p['note']})"); continue
    print(f"{p['phase']:9s} GM {p['GM']:11.2f}  above eq {p['GM_above_equilibrium']:8.2f}"
          f"  driving force {p['driving_force']:9.2f}  state {p['state']}"
          f"  splits? {'splits_if_allowed' in p}")

# --- 7.8 API 3 ------------------------------------------------------------------
r3 = cuni.energy_surface(T=1500, step=0.25)
for p in r3["phases"]:
    print("\n", p["phase"], "approximate:", p["approximate"])
    for xx, g in zip(p["X"], p["GM"]):
        print("   X", xx, " GM", round(g, 1))
r4 = alni.energy_surface(T=1600, step=0.05)
print()
for p in r4["phases"]:
    print(f"{p['phase']:9s} points {p['n_points']:4d}  approximate {p['approximate']}")
