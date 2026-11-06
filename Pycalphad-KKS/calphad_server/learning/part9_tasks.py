# Part 9: end-to-end tasks through the server (start it first:
#   python calphad_server.py --tdb cuni_cost507.tdb --tdb alni_dupin_2001.tdb)
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
from calphad_client import CalphadClient

c = CalphadClient("http://127.0.0.1:8000", db="alni_dupin_2001")

# Task 1: where does Ni-10Al solidify?
r = c.equilibrium({"AL": 0.10}, T=[1700, 1716, 1718, 1719, 1720, 1722], balance="NI")
for res in r["results"]:
    print(res["T"], [(p["phase"], p["state"], p["equivalent_to"], round(p["fraction"], 3),
                      round(p["composition"]["AL"], 4)) for p in res["phases"]],
          "verified:", res["stability_check"]["verified"])

# Task 2: which phases are close to forming at 1718 K?
r2 = c.phase_energies({"AL": 0.10}, T=1718, balance="NI")
print("\nstable:", r2["stable_phases"])
for p in r2["phases"]:
    if p["GM"] is not None:
        print(f"{p['phase']:8s} above eq {p['GM_above_equilibrium']:8.2f}  driving force {p['driving_force']:9.2f}  state {p['state']}")

# Task 3: KKS inputs for FCC_A1 + LIQUID
k = c.kks_parameters(solute="AL", phase_a="FCC_A1", phase_b="LIQUID",
                     T=[1705, 1710, 1715, 1718], mode="parabolic", c_alloy=0.10, Vm=7.0e-6)
print("\nwarnings:", k["warnings"])
for row in k["tie_lines"]:
    print({key: round(row[key], 4) for key in ("T", "c_a", "c_b", "k", "mu", "A_a", "A_b")},
          "consistent:", row["consistent"])
