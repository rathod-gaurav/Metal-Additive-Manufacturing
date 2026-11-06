# Part 6: where pycalphad needs help
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from pycalphad import Database, equilibrium, variables as v
from pycalphad.core.utils import filter_phases

# 1. a missed phase: Ni-20Al-10Cr at 1600 K
db = Database("alcrni.tdb")
comps = ["AL", "CR", "NI", "VA"]
phases = filter_phases(db, comps)
conds = {v.X("AL"): 0.2, v.X("CR"): 0.1, v.T: 1600, v.P: 101325, v.N: 1}
for label, opts in [("default sampling", {}), ("pdens = 2000", {"calc_opts": {"pdens": 2000}})]:
    e = equilibrium(db, comps, phases, conds, **opts).squeeze()
    found = [(str(p), round(float(f), 4)) for p, f in zip(e.Phase.values, e.NP.values) if p != ""]
    print(f"{label:17s}: {found}  GM = {float(e.GM):.2f} J/mol")

# 2. one phase name, two structures: Ni-20Al at 1000 K
db2 = Database("alni_dupin_2001.tdb")
comps2 = ["AL", "NI", "VA"]
e = equilibrium(db2, comps2, filter_phases(db2, comps2), {v.X("AL"): 0.2, v.T: 1000, v.P: 101325, v.N: 1}).squeeze()
for k, p in enumerate(e.Phase.values):
    if p == "":
        continue
    y = e.Y.values[k][:5]       # Y(AL) sub0, Y(NI) sub0, Y(AL) sub1, Y(NI) sub1, Y(VA) sub2
    print(f"{p}: fraction {float(e.NP.values[k]):.4f}, x_Al {float(e.X.values[k][0]):.4f}, "
          f"Al on sublattice 0 = {y[0]:.4f}, Al on sublattice 1 = {y[2]:.4f}")
