# Part 5: pycalphad hands-on
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from pycalphad import Database, equilibrium, calculate, Model, variables as v

db = Database("cuni_cost507.tdb")
comps = ["CU", "NI", "VA"]          # elements + VA (this database uses vacancies)
phases = ["LIQUID", "FCC_A1"]

# ---- 5.3 one equilibrium point
conds = {v.X("NI"): 0.33, v.T: 1500, v.P: 101325, v.N: 1}
eq = equilibrium(db, comps, phases, conds)
print("dimensions:", dict(eq.sizes))
print("variables :", list(eq.data_vars))
e = eq.squeeze()
print("Phase:", e.Phase.values)
print("NP   :", e.NP.values)
print("X    :", e.X.values, "components", list(e.component.values))
print("MU   :", e.MU.values)
print("GM   :", float(e.GM))

# ---- 5.3 many conditions at once
conds = {v.X("NI"): 0.5, v.T: (1500, 1701, 50), v.P: 101325, v.N: 1}
eq = equilibrium(db, comps, phases, conds)
print("Phase array shape:", eq.Phase.shape)
for i, T in enumerate(eq.T.values):
    pt = eq.isel(T=i).squeeze()
    names = [str(p) for p in pt.Phase.values if p != ""]
    fracs = [round(float(f), 4) for f in pt.NP.values if not np.isnan(f)]
    print(f"T = {T:.0f} K: {names} {fracs}")

# ---- 5.4 calculate: G at points you choose
res = calculate(db, comps, "LIQUID", T=1500, P=101325, N=1,
                points=np.array([[0.9, 0.1], [0.5, 0.5], [0.1, 0.9]]), output="GM")
print("calculate dims:", dict(res.sizes))
print("X :", res.X.values.squeeze())
print("GM:", res.GM.values.squeeze())
res = calculate(db, comps, "FCC_A1", T=1500, P=101325, N=1, pdens=50)
print("sampled points for FCC_A1:", res.GM.values.size)

# ---- 5.5 Model: the formula itself
m = Model(db, comps, "LIQUID")
print("site fractions:", m.site_fractions)
print("variables of GM:", sorted(str(s) for s in m.GM.free_symbols))
print("terms:", list(m.models.keys()))
vals = {m.site_fractions[0]: 0.67, m.site_fractions[1]: 0.33, v.T: 1500.0}
print("GM at x_Ni = 0.33, 1500 K:", float(m.GM.xreplace(vals)))
