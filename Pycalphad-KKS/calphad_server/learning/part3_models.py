# Part 3: split G of each phase into its CALPHAD terms (Cu-33Ni, 1500 K)
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
from pycalphad import Database, Model, variables as v

db = Database("cuni_cost507.tdb")
x_ni, T = 0.33, 1500.0
for phase in ("LIQUID", "FCC_A1"):
    m = Model(db, ["CU", "NI", "VA"], phase)
    # site fractions -> numbers (Cu and Ni share the first sublattice; VA fills the second)
    values = {}
    for y in m.site_fractions:
        name = y.species.name
        values[y] = x_ni if name == "NI" else (1 - x_ni if name == "CU" else 1.0)
    values[v.T] = T
    values[v.P] = 101325
    total = 0.0
    print(phase, [str(y) for y in m.site_fractions])
    for key, expr in m.models.items():
        val = float(expr.xreplace(values)) if hasattr(expr, "xreplace") else float(expr)
        total += val
        print(f"   {key:6s} {val:12.1f} J/mol")
    print(f"   total  {total:12.1f}   (m.GM gives {float(m.GM.xreplace(values)):.1f})")
