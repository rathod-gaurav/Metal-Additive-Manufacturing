# Part 4: look inside a TDB database with pycalphad
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
from pycalphad import Database
from tinydb import where

db = Database("cuni_cost507.tdb")
print("elements:", sorted(db.elements))
print("phases:  ", sorted(db.phases))
fcc = db.phases["FCC_A1"]
print("FCC_A1 site ratios:", fcc.sublattices)
print("FCC_A1 constituents:", [sorted(s.name for s in sub) for sub in fcc.constituents])
print("FCC_A1 model hints:", fcc.model_hints)

# every Cu-Ni interaction parameter of the liquid and FCC
query = (where("phase_name").one_of(["LIQUID", "FCC_A1"])) & (where("parameter_type").one_of(["G", "L"]))
for p in db.search(query):
    species = [[s.name for s in sub] for sub in p["constituent_array"]]
    if any(len(sub) > 1 for sub in species):            # more than one species on a sublattice = interaction
        print(p["parameter_type"], p["phase_name"], species, "order", p["parameter_order"], ":", p["parameter"])
