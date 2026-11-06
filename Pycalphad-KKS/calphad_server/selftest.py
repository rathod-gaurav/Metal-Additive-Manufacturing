"""
selftest.py — sanity checks that should hold for ANY database.
    python selftest.py                      # bundled examples
    python selftest.py my.tdb '{"NI":0.1}' NI 1200   # your own: tdb, solutes, balance, T

Checks, at each composition/temperature:
  * phase fractions sum to 1
  * mass balance: sum_p f_p x_p = alloy composition
  * each stable phase lies on the equilibrium tangent plane: G_p = sum_i x_p,i mu_i
  * no single phase is below the equilibrium energy: G_phase >= G_equilibrium
  * stable phases have zero driving force, others <= 0
"""
import json
import sys

import numpy as np

from calphad_core import CalphadSystem

TOL_X, TOL_G = 1e-6, 1e-3  # mole fraction, J/mol


def check(system, solutes, balance, T):
    eq = system.equilibrium(solutes, T, balance=balance)
    r, x = eq["results"], eq["alloy"]
    fr = [p["fraction"] for p in r["phases"]]
    assert abs(sum(fr) - 1) < TOL_X, f"fractions sum to {sum(fr)}"
    for el in x:
        mix = sum(p["fraction"] * p["composition"][el] for p in r["phases"])
        assert abs(mix - x[el]) < TOL_X, f"mass balance {el}: {mix} vs {x[el]}"
    for p in r["phases"]:
        plane = sum(p["composition"][el] * r["chemical_potentials"][el] for el in x)
        assert abs(plane - p["GM"]) < TOL_G
    g_sys = sum(p["fraction"] * p["GM"] for p in r["phases"])
    assert abs(g_sys - r["GM_system"]) < TOL_G

    pe = system.phase_energies(solutes, T, balance=balance)
    stable = {p["phase"] for p in r["phases"]}
    for p in pe["phases"]:
        if p["GM"] is None:
            continue
        assert p["GM"] >= pe["GM_equilibrium"] - TOL_G, f"{p['phase']} below equilibrium"
        df = p.get("driving_force")
        if df is not None:
            assert df <= 1e-3, f"{p['phase']} positive driving force {df}"
            if p["phase"] in stable and len(stable) == len(r["phases"]):
                assert abs(df) < 1e-2, f"stable {p['phase']} has driving force {df}"
    names = ", ".join(f"{p['phase']} {p['fraction']:.3f}" for p in r["phases"])
    print(f"  ok  {system.name:14s} T={T:6.0f}  {json.dumps(solutes)} bal {balance}: {names}")


if __name__ == "__main__":
    if len(sys.argv) == 5:
        cases = [(sys.argv[1], json.loads(sys.argv[2]), sys.argv[3], float(sys.argv[4]))]
    else:
        cases = [("cuni_cost507.tdb", {"NI": x}, "CU", T)
                 for x in (0.05, 0.33, 0.6, 0.95) for T in (1400, 1500, 1600, 1700)]
        cases += [("alcrni.tdb", {"AL": a, "CR": c}, "NI", T)
                  for a, c in ((0.1, 0.1), (0.2, 0.1), (0.25, 0.1), (0.4, 0.2))
                  for T in (900, 1273, 1600)]
    systems = {}
    for tdb, sol, bal, T in cases:
        if tdb not in systems:
            systems[tdb] = CalphadSystem.from_file(tdb)
        check(systems[tdb], sol, bal, T)
    print(f"all {len(cases)} cases passed")
