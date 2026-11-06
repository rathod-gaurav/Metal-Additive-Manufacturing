# Part 7: calling the running server from Python (start the server first)
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import requests

URL = "http://127.0.0.1:8000"
body = {"db": "cuni_cost507", "composition": {"NI": 0.33}, "balance": "CU", "T": 1500}
resp = requests.post(URL + "/equilibrium", json=body)
print("status:", resp.status_code)
data = resp.json()                      # JSON text -> Python dicts and lists
for ph in data["results"]["phases"]:
    print(ph["phase"], round(ph["fraction"], 4), round(ph["composition"]["NI"], 4))

# the same through the ready-made client (standard library only)
from calphad_client import CalphadClient
c = CalphadClient(URL, db="cuni_cost507")
r = c.equilibrium({"NI": 0.33}, T=1500, balance="CU")
print("client:", [(p["phase"], round(p["fraction"], 4)) for p in r["results"]["phases"]])
try:
    c.equilibrium({"NI": 0.3, "FE": 0.7}, T=1500)
except Exception as exc:
    print(type(exc).__name__, "->", exc)
