# Part 2: common tangent, chemical potentials, tie line and lever rule for Cu-33Ni at 1500 K
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from calphad_core import CalphadSystem

s = CalphadSystem.from_file("cuni_cost507.tdb")
r = s.equilibrium({"NI": 0.33}, T=1500, balance="CU")["results"]
mu = r["chemical_potentials"]
print("chemical potentials:", {k: round(v, 2) for k, v in mu.items()})
for p in r["phases"]:
    print(p["phase"], "fraction", round(p["fraction"], 4), "x_Ni", round(p["composition"]["NI"], 4))

# the common tangent line: G_line(x) = mu_Cu * (1 - x) + mu_Ni * x
line = lambda x: mu["CU"] * (1 - x) + mu["NI"] * x
print("slope of the tangent (mu_Ni - mu_Cu):", round(mu["NI"] - mu["CU"], 2), "J/mol")

# lever rule by hand
x0 = 0.33
xL = next(p["composition"]["NI"] for p in r["phases"] if p["phase"] == "LIQUID")
xS = next(p["composition"]["NI"] for p in r["phases"] if p["phase"] == "FCC_A1")
f_liquid = (xS - x0) / (xS - xL)
print(f"lever rule: liquid fraction = ({xS:.4f} - {x0}) / ({xS:.4f} - {xL:.4f}) = {f_liquid:.4f}")

# G of each phase minus the tangent line, for the chart
surf = s.energy_surface(T=1500, step=0.02)
i_ni = surf["elements"].index("NI")
for ph in surf["phases"]:
    X = np.array(ph["X"])[:, i_ni]; G = np.array(ph["GM"])
    keep = (X >= 0.18 - 1e-9) & (X <= 0.50 + 1e-9)
    print(ph["phase"], [(round(x, 2), round(g - line(x), 1)) for x, g in zip(X[keep], G[keep])])
