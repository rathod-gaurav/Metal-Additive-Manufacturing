# Part 7: pycalphad's Workspace, IsolatedPhase and DormantPhase (what API 1 and API 2 use)
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from pycalphad import Database, Workspace, variables as v
from pycalphad.property_framework.metaproperties import IsolatedPhase, DormantPhase

db = Database("alni_dupin_2001.tdb")
comps = ["AL", "NI", "VA"]
phases = ["FCC_L12", "BCC_B2", "LIQUID"]
wks = Workspace(db, comps, phases, {v.X("AL"): 0.25, v.T: 1600, v.P: 101325, v.N: 1})

print("GM of the equilibrium:", float(np.ravel(wks.get("GM"))[0]))
print("mu_Al, mu_Ni:", [round(float(np.ravel(wks.get(q))[0]), 1) for q in ("MU(AL)", "MU(NI)")])
liq = IsolatedPhase("LIQUID", wks=wks)
print("G of LIQUID alone at this composition:", round(float(np.ravel(wks.get(liq("GM")))[0]), 1))
print("driving force of LIQUID:", round(float(np.ravel(wks.get(DormantPhase("LIQUID", wks=wks).driving_force))[0]), 1))
