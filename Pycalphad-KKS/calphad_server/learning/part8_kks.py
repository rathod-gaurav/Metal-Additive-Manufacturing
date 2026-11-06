# Part 8: KKS building blocks (run from inside the calphad_server folder)
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)   # find the framework and the .tdb files
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from calphad_kks import KKSModel

# --- 8.2 interface parameters ------------------------------------------------
eps2, W = KKSModel.interface_parameters(sigma=0.3, width=25e-9)
print(f"eps^2 = {eps2:.3e} J/m,  W = {W:.3e} J/m^3")
print(f"check sigma = {np.sqrt(eps2)*np.sqrt(W)/(3*np.sqrt(2)):.4f} J/m^2,"
      f" width = {2.2*np.sqrt(2)*np.sqrt(eps2)/np.sqrt(W)*1e9:.2f} nm")

# --- 8.5 the KKS model for Cu-Ni, solid = FCC_A1 at phi = 1 --------------------
kks = KKSModel.from_tdb("cuni_cost507.tdb", solute="NI", phase_a="FCC_A1",
                        phase_b="LIQUID", mode="exact", Vm=7.1e-6)
row = kks.tie_line(1500)
print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})

# --- 8.3 naive mixture vs KKS in the middle of an equilibrium interface ---------
cS, cL, mu = row["c_a"], row["c_b"], row["mu"]
h = 0.5
c = h * cS + (1 - h) * cL
tangent = lambda cc: row["G_b"] + mu * (cc - cL)          # the common tangent line
(GS, _, _), (GL, _, _) = kks.energies(c, c, 1500)        # both phases at the same c
naive = h * GS + (1 - h) * GL - tangent(c)
(GS2, _, _), (GL2, _, _) = kks.energies(cS, cL, 1500)    # KKS: each at its own c
kks_e = h * GS2 + (1 - h) * GL2 - tangent(c)
print(f"mid-interface c = {float(c):.4f}: energy above the tangent, naive {float(naive):.2f} J/mol,"
      f" KKS {float(kks_e):.2e} J/mol")

# --- 8.8 a static interface -------------------------------------------------------
x = np.linspace(-3, 3, 7)                      # position, in units of the interface width
phi = 0.5 * (1 - np.tanh(x))                   # solid on the left, liquid on the right
hh = kks.h(phi)
c = hh * cS + (1 - hh) * cL                    # equilibrium composition profile
for T in (1500, 1490, 1510):
    ca, cb, mu = kks.split(c, phi=phi, T=T)
    dF = kks.driving_force(ca, cb, mu, T)
    print(f"\nT = {T} K")
    print("   phi      c     c_S    c_L    mu (J/mol)   dF (J/m^3)")
    for i in range(len(x)):
        print(f"{phi[i]:6.3f} {c[i]:7.4f} {ca[i]:6.4f} {cb[i]:6.4f} {mu[i]:11.1f} {dF[i]:12.4g}")
