"""
kks_1d_directional.py
=====================
1-D directional solidification of Ni-10Al (at.%) with the KKS model coupled to
CALPHAD (Dupin 2001 Al-Ni database), in a moving temperature gradient

    T(x, t) = T_sol + G (x - x0 - V t)          (frozen-temperature approximation)

the setting used for melt pools in additive manufacturing. Shows how
calphad_kks plugs into a solver loop with a temperature FIELD.

Physics checks (planar front, steady state, local equilibrium):
  * solute is conserved to round-off
  * the solid forms at the alloy composition c0
  * the interface settles just below the solidus temperature of c0
    (below by the kinetic undercooling of the finite interface mobility)
  * the liquid at the interface reaches about c0 / k

Run:  python examples/kks_1d_directional.py [--mode exact|parabolic|table]
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from calphad_kks import KKSModel  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--mode", default="exact", choices=["exact", "parabolic", "table"])
ap.add_argument("--tdb", default=os.path.join(os.path.dirname(HERE), "alni_dupin_2001.tdb"))
ap.add_argument("--out", default=HERE)
args = ap.parse_args()

# ---- material and process parameters (illustrative, SI units) -------------
c0 = 0.10          # Ni-10Al, mole fraction Al
Vm = 7.0e-6        # m^3/mol
sigma = 0.3        # J/m^2, solid-liquid interface energy (illustrative)
DL, DS = 3e-9, 1e-13   # m^2/s
G = 1.0e7          # K/m, thermal gradient (AM melt pools: ~1e6-1e7)
V = 0.01           # m/s, isotherm velocity
N, dx = 600, 5e-9  # 3 um domain
width = 5 * dx     # interface width 2*lambda
M_factor = 1.0     # phase-field mobility, numerical choice (see calphad_kks docs)
t_end = 1.8e-4     # s

# ---- thermodynamics from CALPHAD -------------------------------------------
temps = np.arange(1690.0, 1725.0 + 1e-9, 1.0)
kks = KKSModel.from_tdb(args.tdb, solute="AL", phase_a="FCC_A1", phase_b="LIQUID",
                        mode=args.mode, Vm=Vm, c_alloy=c0,
                        T=None if args.mode == "exact" else temps)
for w in kks.warnings:
    print("warning:", w)
ref = KKSModel.from_tdb(args.tdb, solute="AL", phase_a="FCC_A1", phase_b="LIQUID",
                        mode="exact", c_alloy=c0)
rows = [r for r in ref.tie_lines(np.arange(1700.0, 1722.0, 0.25)) if r]
Tg = np.array([r["T"] for r in rows])
ca_eq = np.array([r["c_a"] for r in rows])
cb_eq = np.array([r["c_b"] for r in rows])
T_sol = float(np.interp(c0, ca_eq[::-1], Tg[::-1]))       # solidus of c0
T_liq = float(np.interp(c0, cb_eq[::-1], Tg[::-1]))       # liquidus of c0
k_sol = float(np.interp(T_sol, Tg, ca_eq / cb_eq))
print(f"Ni-{100*c0:.0f}Al: solidus {T_sol:.2f} K, liquidus {T_liq:.2f} K, "
      f"k(T_sol) = {k_sol:.4f}, c0/k = {c0/k_sol:.4f}")

# ---- phase-field parameters --------------------------------------------------
eps2, W = KKSModel.interface_parameters(sigma, width)
M = M_factor * DL / eps2
dt = 0.2 * dx ** 2 / DL
nsteps = int(t_end / dt)
x = (np.arange(N) + 0.5) * dx
x0 = 0.4e-6

phi = 0.5 * (1 - np.tanh(np.sqrt(2 * W) * (x - x0) / (2 * np.sqrt(eps2))))
c = np.full(N, c0)
mass0 = c.sum() * dx
ca, cb = c.copy(), c.copy()

hist_t, hist_xi, hist_Ti, hist_cs, hist_cl = [], [], [], [], []
snaps = []
t0 = time.time()
for n in range(nsteps + 1):
    t = n * dt
    T = T_sol + G * (x - x0 - V * t)                     # temperature field
    h = kks.h(phi)
    ca, cb, mu = kks.split(c, phi=phi, T=T, ca0=ca, cb0=cb)   # warm start
    dF = kks.driving_force(ca, cb, mu, T)                # J/m^3

    if n % 2000 == 0 or n == nsteps:
        i = np.flatnonzero((phi[:-1] >= 0.5) & (phi[1:] < 0.5))
        if len(i):
            i = i[-1]
            f = (phi[i] - 0.5) / (phi[i] - phi[i + 1])
            xi = x[i] + f * dx
            hist_t.append(t); hist_xi.append(xi)
            hist_Ti.append(T_sol + G * (xi - x0 - V * t))
            hist_cs.append(np.interp(xi - 0.3e-6, x, c))
            hist_cl.append(cb[i + 1])
    if n in (0, nsteps // 3, 2 * nsteps // 3, nsteps):
        snaps.append((t, phi.copy(), c.copy(), T.copy()))
    if n == nsteps:
        break

    lap = np.empty(N)
    lap[1:-1] = phi[2:] - 2 * phi[1:-1] + phi[:-2]
    lap[0], lap[-1] = phi[1] - phi[0], phi[-2] - phi[-1]
    dphi = M * (eps2 * lap / dx ** 2 - W * kks.dg(phi) - kks.dh(phi) * dF)

    hf = 0.5 * (h[1:] + h[:-1])
    Df = hf * DS + (1 - hf) * DL
    flux = -Df * (hf * np.diff(ca) + (1 - hf) * np.diff(cb)) / dx
    dc = np.zeros(N)
    dc[:-1] -= flux / dx
    dc[1:] += flux / dx

    phi = np.clip(phi + dt * dphi, 0.0, 1.0)
    c = c + dt * dc

wall = time.time() - t0
drift = (c.sum() * dx - mass0) / mass0
print(f"mode={args.mode}: {nsteps} steps in {wall:.1f} s ({1e6*wall/nsteps:.0f} us/step)")
print(f"solute drift {drift:.1e}")
print(f"final: interface T {hist_Ti[-1]:.2f} K (solidus {T_sol:.2f}), "
      f"solid c {hist_cs[-1]:.4f} (c0 {c0}), interface liquid c {hist_cl[-1]:.4f} "
      f"(c0/k {c0/k_sol:.4f}), interface speed "
      f"{(hist_xi[-1]-hist_xi[-4])/(hist_t[-1]-hist_t[-4]):.4f} m/s (V {V})")

np.savez(os.path.join(args.out, f"kks_1d_{args.mode}.npz"), t=hist_t, xi=hist_xi,
         Ti=hist_Ti, cs=hist_cs, cl=hist_cl, T_sol=T_sol, T_liq=T_liq, c0=c0, k=k_sol)

# ---- figure ------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#b9d3f2", "#7aa9e3", "#3f86d9", "#143c6b"]


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), facecolor="#fcfcfb")
ax = axs[0]; style(ax)
for (t, ph, cc, TT), col in zip(snaps, RAMP):
    ax.plot(x * 1e6, cc, color=col, lw=1.8, label=f"t = {t*1e6:.0f} µs")
ax.axhline(c0, color=INK2, lw=1, ls="--")
ax.axhline(c0 / k_sol, color=ORANGE, lw=1, ls=":")
ax.text(0.05, c0, " c0", color=INK2, fontsize=8, va="bottom")
ax.text(x[-1] * 1e6, c0 / k_sol, "c0/k ", color=INK2, fontsize=8, va="bottom", ha="right")
ax.set_xlabel("position  [µm]", color=INK); ax.set_ylabel("x_Al", color=INK)
ax.set_title("Composition (front moves right)", loc="left", fontsize=10, color=INK)
ax.legend(fontsize=7.5, frameon=False, loc="upper left")

ax = axs[1]; style(ax)
ax.plot(np.array(hist_t) * 1e6, hist_Ti, color=BLUE, lw=2, label="interface T")
ax.axhline(T_sol, color=INK2, lw=1, ls="--")
ax.axhline(T_liq, color=INK2, lw=1, ls=":")
ax.text(1, T_sol, "solidus of c0", fontsize=8, color=INK2, va="bottom")
ax.text(1, T_liq, "liquidus of c0", fontsize=8, color=INK2, va="bottom")
ax.set_xlabel("time  [µs]", color=INK); ax.set_ylabel("T  [K]", color=INK)
ax.set_title("Interface temperature → solidus (steady state)", loc="left", fontsize=10, color=INK)

ax = axs[2]; style(ax)
ax.plot(np.array(hist_t) * 1e6, hist_cs, color=BLUE, lw=2, label="solid, 0.3 µm behind front")
ax.plot(np.array(hist_t) * 1e6, hist_cl, color=ORANGE, lw=2, label="liquid at the front (c_L)")
ax.axhline(c0, color=INK2, lw=1, ls="--")
ax.axhline(c0 / k_sol, color=INK2, lw=1, ls=":")
ax.set_xlabel("time  [µs]", color=INK); ax.set_ylabel("x_Al", color=INK)
ax.set_title("Solid → c0, interface liquid → c0/k", loc="left", fontsize=10, color=INK)
ax.legend(fontsize=7.5, frameon=False, loc="center right")
fig.suptitle(f"Ni–10Al directional solidification, KKS + CALPHAD ({args.mode} energies), "
             f"G = {G:.0e} K/m, V = {V} m/s", x=0.01, ha="left", fontsize=11, color=INK)
fig.tight_layout()
path = os.path.join(args.out, f"kks_1d_{args.mode}.png")
fig.savefig(path, dpi=150, facecolor=fig.get_facecolor())
print("figure:", path)
