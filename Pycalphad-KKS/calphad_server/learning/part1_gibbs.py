# Part 1: Gibbs energy of pure Cu (solid vs liquid) and the ideal entropy of mixing
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from pycalphad import Database, calculate

db = Database("cuni_cost507.tdb")
temps = [1300, 1358, 1400, 1500]
for T in temps:
    g_fcc = calculate(db, ["CU", "NI", "VA"], "FCC_A1", T=T, P=101325, N=1,
                      points=np.array([[1 - 1e-12, 1e-12, 1.0]]), output="GM").GM.values.item()
    g_liq = calculate(db, ["CU", "NI", "VA"], "LIQUID", T=T, P=101325, N=1,
                      points=np.array([[1 - 1e-12, 1e-12]]), output="GM").GM.values.item()
    print(f"T = {T} K: G_fcc = {g_fcc:10.1f} | G_liquid = {g_liq:10.1f} |  G_liquid - G_fcc = {g_liq - g_fcc:7.1f} J/mol")

# where do they cross? (melting point of pure Cu in this database)
Ts = np.arange(1300, 1420, 0.5)
pts_f = np.array([[1 - 1e-12, 1e-12, 1.0]]); pts_l = np.array([[1 - 1e-12, 1e-12]])
gf = np.array([calculate(db, ["CU", "NI", "VA"], "FCC_A1", T=t, P=101325, N=1, points=pts_f, output="GM").GM.values.item() for t in Ts])
gl = np.array([calculate(db, ["CU", "NI", "VA"], "LIQUID", T=t, P=101325, N=1, points=pts_l, output="GM").GM.values.item() for t in Ts])
i = np.argmin(abs(gl - gf)); print("crossing near", Ts[i], "K")

R = 8.314462618
def g_mix_ideal(x, T):
    """Ideal Gibbs energy of mixing, J/mol: R T [x ln x + (1-x) ln(1-x)]"""
    return R * T * (x * np.log(x) + (1 - x) * np.log(1 - x))

for x in (1e-12, 0.1, 0.25, 0.5):
    print(f"ideal mixing at x = {x}, 1500 K: {g_mix_ideal(x, 1500):.1f} J/mol")
