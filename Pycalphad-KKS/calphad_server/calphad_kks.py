"""
calphad_kks.py
==============
Kim-Kim-Suzuki (KKS) building blocks for a binary A-B alloy and two phases,
from any TDB database. Sits on top of calphad_core (which it does not change).

Inside a phase-field code (NumPy arrays of any shape; T may be a field):

    from calphad_kks import KKSModel
    kks = KKSModel.from_tdb("alni_dupin_2001.tdb", solute="AL",
                            phase_a="FCC_A1", phase_b="LIQUID",   # phi = 1 is phase_a
                            mode="exact", Vm=7.0e-6, c_alloy=0.10)
    ca, cb, mu = kks.split(c, phi, T)          # KKS constraint at every grid point
    dF = kks.driving_force(ca, cb, mu, T)      # J/m^3, the term in the phi equation
    # dphi/dt = M [ eps^2 lap(phi) - W g'(phi) - h'(phi) dF ]
    # dc/dt   = div[ D(phi) ( h grad(ca) + (1-h) grad(cb) ) ]

Free-energy representations (`mode`)
-------------------------------------
  "exact"     The CALPHAD expression compiled to machine code, with exact first
              and second derivatives. Substitutional phases only (one mixing
              sublattice, e.g. LIQUID, FCC_A1, BCC_A2 without vacancies).
  "table"     G(c, T) computed from CALPHAD on a grid and splined (the ideal
              mixing term is handled analytically). Works for any phase,
              including ordered ones (each grid point is the phase minimised on
              its own at that composition).
  "parabolic" f = G_eq + mu_eq (c - c_eq) + A/2 (c - c_eq)^2 about the tie line
              at each temperature (A = curvature from CALPHAD). The usual KKS
              simplification; exact at the tie line, approximate away from it.

Conventions: c = mole fraction of the solute B. Energies in J/mol of atoms
unless a method says per volume (J/m^3 = J/mol / Vm). mu is the diffusion
potential dG/dc = mu_B - mu_A.
"""

import json

import numpy as np
import symengine as se
from pycalphad import Workspace, equilibrium, variables as v
from scipy.interpolate import CubicSpline, RectBivariateSpline

from calphad_core import P_DEFAULT, CalphadError, CalphadSystem

R = 8.314462618  # J/(mol K)
EPS = 1e-12


def _bcast(c, T):
    c = np.asarray(c, dtype=float)
    T = np.asarray(T, dtype=float)
    c, T = np.broadcast_arrays(c, T)
    return c, T


def _ideal(c, T):
    """R T [c ln c + (1-c) ln(1-c)] and its first two c-derivatives."""
    c = np.clip(c, EPS, 1 - EPS)
    g = R * T * (c * np.log(c) + (1 - c) * np.log(1 - c))
    d1 = R * T * (np.log(c) - np.log(1 - c))
    d2 = R * T / (c * (1 - c))
    return g, d1, d2


# ===========================================================================
# Phase free-energy representations: __call__(c, T) -> (G, dG/dc, d2G/dc2)
# ===========================================================================
class PhaseEnergy:
    kind = "base"

    def __init__(self, name):
        self.name = name

    def __call__(self, c, T):
        raise NotImplementedError


class ExactPhase(PhaseEnergy):
    """CALPHAD G(c, T) compiled from the database. Substitutional phases only."""
    kind = "exact"

    def __init__(self, system, phase, solute, solvent):
        super().__init__(phase)
        comps = system._components([solute, solvent])
        model = system._model(comps, phase)
        if system._mixing_sublattice(model, [solute, solvent]) is None:
            raise CalphadError(
                f"{phase} has sublattices or vacancies mixing with the elements, so G at "
                "fixed composition needs an internal minimisation. Use mode='table'.")
        x = se.Symbol("x")
        subs = {}
        for y in model.site_fractions:
            name = y.species.name
            subs[y] = x if name == solute else (1 - x if name == solvent else 1)
        G = model.GM.xreplace(subs).xreplace({v.P: P_DEFAULT})
        dG = G.diff(x)
        d2G = dG.diff(x)
        backend = "llvm" if getattr(se, "have_llvm", False) else "lambda"
        self._fn = se.Lambdify([x, v.T], [G, dG, d2G], backend=backend, cse=True)

    def __call__(self, c, T):
        c, T = _bcast(c, T)
        shape = c.shape
        cc = np.clip(c.ravel(), EPS, 1 - EPS)
        out = np.asarray(self._fn(np.stack([cc, T.ravel()], axis=-1))).reshape(cc.size, 3)
        return (out[:, 0].reshape(shape), out[:, 1].reshape(shape), out[:, 2].reshape(shape))


class TabulatedPhase(PhaseEnergy):
    """G(c, T) from CALPHAD on a grid, splined. Any phase.

    The ideal-mixing term R T [c ln c + (1-c) ln(1-c)] is subtracted before
    splining and added back analytically, so the log singularities at the
    composition ends do not spoil the fit.
    """
    kind = "table"

    def __init__(self, system, phase, solute, solvent, T_grid, c_grid=None,
                 n_c=201, c_range=(1e-4, 1 - 1e-4), branch="lowest"):
        """branch (order/disorder phases only):
             "lowest"     the lower of the ordered and disordered states at each
                          point (G has a kink where they swap)
             "ordered"    the ordered branch only (e.g. gamma-prime), including
                          where it is metastable; smooth, which suits KKS
             "disordered" the disordered branch only
        """
        super().__init__(phase)
        if branch not in ("lowest", "ordered", "disordered"):
            raise CalphadError("branch must be 'lowest', 'ordered' or 'disordered'")
        T_grid = np.atleast_1d(np.asarray(T_grid, dtype=float))
        c_grid = (np.linspace(c_range[0], c_range[1], n_c) if c_grid is None
                  else np.asarray(c_grid, dtype=float))
        comps = system._components([solute, solvent])
        model = system._model(comps, phase)
        ordered = np.zeros((len(T_grid), len(c_grid)), dtype=bool)
        with system._lock:
            if system._mixing_sublattice(model, [solute, solvent]) is not None:
                ex = ExactPhase(system, phase, solute, solvent)
                CC, TT = np.meshgrid(c_grid, T_grid)
                G = ex(CC, TT)[0]
            else:
                # Minimise the phase on its own at every grid point, from several
                # starting configurations (disordered + ordered), so both branches
                # of an order/disorder phase are found. Vectorised over T.
                nT, nc = len(T_grid), len(c_grid)
                starts0 = system._starting_points(model, {solvent: 0.5, solute: 0.5})
                G_all = np.full((len(starts0), nT, nc), np.nan)
                ord_all = np.zeros((len(starts0), nT, nc), dtype=bool)
                sfs = model.site_fractions
                subl = sorted({sf.sublattice_index for sf in sfs if sf.species.name == solute})
                for j, cj in enumerate(c_grid):
                    x = {solvent: 1.0 - cj, solute: cj}
                    for k, y0 in enumerate(system._starting_points(model, x)):
                        try:
                            GM, _, Y, _ = system._refine(comps, phase, x, T_grid, P_DEFAULT, y0)
                        except Exception:
                            continue
                        G_all[k, :, j] = GM
                        if len(subl) > 1:
                            ys = [Y[:, [i for i, sf in enumerate(sfs)
                                        if sf.sublattice_index == s_ and sf.species.name == solute][0]]
                                  for s_ in subl]
                            spread = np.max([np.abs(yy - ys[0]) for yy in ys[1:]], axis=0)
                            ord_all[k, :, j] = np.nan_to_num(spread) > 1e-3
                with np.errstate(all="ignore"):
                    G_ord = np.nanmin(np.where(ord_all, G_all, np.nan), axis=0)
                    G_dis = np.nanmin(np.where(~ord_all & np.isfinite(G_all), G_all, np.nan), axis=0)
                    G_low = np.fmin(G_ord, G_dis)
                ordered = np.isfinite(G_ord) & ~(G_dis < G_ord)
                G = {"lowest": G_low, "ordered": G_ord, "disordered": G_dis}[branch]
                if branch == "ordered":
                    ordered = np.isfinite(G_ord)
        ok = np.all(np.isfinite(G), axis=0)          # compositions valid at every T
        if ok.sum() < 8:
            raise CalphadError(f"{phase}: too few usable points to tabulate "
                               f"({ok.sum()} of {len(c_grid)}, branch={branch}).")
        # keep the longest contiguous run of usable compositions
        runs, start = [], None
        for i, flag in enumerate(np.append(ok, False)):
            if flag and start is None:
                start = i
            elif not flag and start is not None:
                runs.append((start, i - 1))
                start = None
        first, last = max(runs, key=lambda r: r[1] - r[0])
        if last - first + 1 < 8:
            raise CalphadError(f"{phase}: usable composition range too short to tabulate.")
        sub = ordered[:, first:last + 1]
        self.branch = branch
        self.mixed_branches = bool(branch == "lowest" and sub.any() and (~sub).any())
        self.c_grid = c_grid[first:last + 1]
        self.T_grid = T_grid
        self.G_table = G[:, first:last + 1]
        CC, TT = np.meshgrid(self.c_grid, T_grid)
        Gex = self.G_table - _ideal(CC, TT)[0]
        if len(T_grid) == 1:
            self._spl1 = CubicSpline(self.c_grid, Gex[0])
            self._spl2 = None
        else:
            self._spl1 = None
            self._spl2 = RectBivariateSpline(self.c_grid, T_grid, Gex.T,
                                             kx=3, ky=min(3, len(T_grid) - 1))
        self.c_min, self.c_max = float(self.c_grid[0]), float(self.c_grid[-1])

    def __call__(self, c, T):
        c, T = _bcast(c, T)
        cc = np.clip(c, self.c_min, self.c_max)
        if self._spl1 is not None:
            g0 = self._spl1(cc)
            g1 = self._spl1(cc, 1)
            g2 = self._spl1(cc, 2)
        else:
            Tc = np.clip(T, self.T_grid[0], self.T_grid[-1])
            g0 = self._spl2.ev(cc, Tc)
            g1 = self._spl2.ev(cc, Tc, dx=1)
            g2 = self._spl2.ev(cc, Tc, dx=2)
        i0, i1, i2 = _ideal(cc, T)
        return g0 + i0, g1 + i1, g2 + i2


class ParabolicPhase(PhaseEnergy):
    """f = G_eq + mu_eq (c - c_eq) + A/2 (c - c_eq)^2 with coefficients vs T.

    Both phases of a KKSModel share the tie line's mu_eq and touch the same
    common tangent, so equilibrium is reproduced exactly at every tabulated T
    (and to interpolation accuracy in between).
    """
    kind = "parabolic"

    def __init__(self, name, T_table, c_eq, A, mu_eq, G_eq):
        super().__init__(name)
        self.T_table = np.atleast_1d(np.asarray(T_table, dtype=float))
        self.c_eq = np.atleast_1d(np.asarray(c_eq, dtype=float))
        self.A = np.atleast_1d(np.asarray(A, dtype=float))
        self.mu_eq = np.atleast_1d(np.asarray(mu_eq, dtype=float))
        self.G_eq = np.atleast_1d(np.asarray(G_eq, dtype=float))

    def coefficients(self, T):
        T = np.asarray(T, dtype=float)
        if len(self.T_table) == 1:
            f = lambda arr: np.full_like(T, arr[0])  # noqa: E731
        else:
            f = lambda arr: np.interp(T, self.T_table, arr)  # noqa: E731
        return f(self.c_eq), f(self.A), f(self.mu_eq), f(self.G_eq)

    def __call__(self, c, T):
        c, T = _bcast(c, T)
        ce, A, mu, G0 = self.coefficients(T)
        d = c - ce
        return G0 + mu * d + 0.5 * A * d * d, mu + A * d, A


# ===========================================================================
# The KKS model
# ===========================================================================
class KKSModel:
    """Two-phase KKS model for a binary alloy.

    Parameters
    ----------
    system   : CalphadSystem (or use KKSModel.from_tdb)
    solute   : element B; c = x_B
    phase_a  : phase at phi = 1 (e.g. the solid, "FCC_A1")
    phase_b  : phase at phi = 0 (e.g. "LIQUID")
    mode     : "exact" | "table" | "parabolic"
    Vm       : molar volume, m^3/mol (needed for anything per volume)
    c_alloy  : nominal alloy composition. Picks the right two-phase field when
               the same two phases meet in more than one place.
    T        : temperatures for "table"/"parabolic" (list), or (Tmin, Tmax, step).
               Not needed for "exact".
    solvent  : element A (defaults to the other element of a binary database)
    reference: "exact" | "table" — what the parabolic fit and the tie lines
               are computed from (default: exact when both phases allow it)
    branches : for order/disorder phases in "table"/"parabolic" mode, e.g.
               {"FCC_L12": "ordered"} for gamma-prime (see TabulatedPhase)

    Check `warnings` after construction, and the "consistent" flag of each tie
    line (False = the KKS energies disagree with pycalphad's own equilibrium).
    """

    def __init__(self, system, solute, phase_a, phase_b, mode="exact", Vm=None,
                 c_alloy=None, T=None, solvent=None, reference=None, n_c=201,
                 c_range=(1e-4, 1 - 1e-4), branches=None):
        self.system = system
        self.B = solute.upper()
        others = [e for e in system.elements if e != self.B]
        if solvent is None:
            if len(others) != 1:
                raise CalphadError(f"Database has {system.elements}; name the solvent.")
            solvent = others[0]
        self.A = solvent.upper()
        system._check_elements([self.A, self.B])
        self.phase_names = (phase_a.upper(), phase_b.upper())
        for ph in self.phase_names:
            if ph not in system.db.phases:
                raise CalphadError(f"Unknown phase {ph}. Database has {system.phases}.")
        if mode not in ("exact", "table", "parabolic"):
            raise CalphadError("mode must be 'exact', 'table' or 'parabolic'")
        self.mode = mode
        self.Vm = Vm
        self.c_alloy = c_alloy
        self.comps = system._components([self.A, self.B])
        self._tie_cache = {}

        T_list = self._temperature_list(T)
        if mode in ("table", "parabolic") and T_list is None:
            raise CalphadError(f"mode='{mode}' needs temperatures: T=[...] or (Tmin, Tmax, step).")

        # reference energies: exact where possible, otherwise tabulated
        branches = {k.upper(): val for k, val in (branches or {}).items()}
        ref = []
        for ph in self.phase_names:
            want = reference or ("table" if mode == "table" else "exact")
            if want == "exact":
                try:
                    ref.append(ExactPhase(system, ph, self.B, self.A))
                    continue
                except CalphadError:
                    if mode == "exact":
                        raise
                    if T_list is None:
                        raise
            ref.append(TabulatedPhase(system, ph, self.B, self.A, T_list, n_c=n_c,
                                      c_range=c_range, branch=branches.get(ph, "lowest")))
        self.reference = tuple(ref)
        self.warnings = []
        for r in ref:
            if getattr(r, "mixed_branches", False):
                self.warnings.append(
                    f"{r.name} switches between ordered and disordered states inside the "
                    "tabulated range, so G has a kink there. Pass branches={'"
                    f"{r.name}': 'ordered'}} (or use the disordered phase, e.g. FCC_A1).")

        if mode == "parabolic":
            rows = self.tie_lines(T_list)
            good = [r for r in rows if r is not None]
            if not good:
                raise CalphadError(f"{self.phase_names} never coexist in {T_list[0]}-{T_list[-1]} K.")
            Ts = [r["T"] for r in good]
            pa = ParabolicPhase(self.phase_names[0], Ts, [r["c_a"] for r in good],
                                [r["A_a"] for r in good], [r["mu"] for r in good],
                                [r["G_a"] for r in good])
            pb = ParabolicPhase(self.phase_names[1], Ts, [r["c_b"] for r in good],
                                [r["A_b"] for r in good], [r["mu"] for r in good],
                                [r["G_b"] for r in good])
            self.phases = (pa, pb)
            self.T_valid = (min(Ts), max(Ts))
        else:
            self.phases = self.reference
            self.T_valid = None if T_list is None else (float(T_list[0]), float(T_list[-1]))

    @classmethod
    def from_tdb(cls, path, solute, phase_a, phase_b, **kw):
        return cls(CalphadSystem.from_file(path), solute, phase_a, phase_b, **kw)

    @staticmethod
    def _temperature_list(T):
        if T is None:
            return None
        if isinstance(T, tuple) and len(T) == 3:
            lo, hi, step = T
            return np.arange(lo, hi + 0.5 * step, step, dtype=float)
        return np.sort(np.atleast_1d(np.asarray(T, dtype=float)))

    # -- interpolation functions --------------------------------------------------
    @staticmethod
    def h(phi):
        return phi ** 3 * (6 * phi ** 2 - 15 * phi + 10)

    @staticmethod
    def dh(phi):
        return 30 * phi ** 2 * (1 - phi) ** 2

    @staticmethod
    def g(phi):
        return phi ** 2 * (1 - phi) ** 2

    @staticmethod
    def dg(phi):
        return 2 * phi * (1 - phi) * (1 - 2 * phi)

    @staticmethod
    def interface_parameters(sigma, width):
        """eps^2 [J/m] and W [J/m^3] from interface energy sigma [J/m^2] and
        width 2*lambda [m] (phi from 0.1 to 0.9), for
        f = ... + W g(phi) + (eps^2/2)|grad phi|^2 (KKS 1999):
            sigma = eps sqrt(W) / (3 sqrt 2),  2 lambda = 2.2 sqrt(2) eps / sqrt(W)."""
        eps2 = 3.0 * sigma * width / 2.2
        W = 13.2 * sigma / width
        return eps2, W

    # -- phase energies ---------------------------------------------------------------
    def energies(self, ca, cb, T):
        """(G_a, dG_a, d2G_a), (G_b, dG_b, d2G_b) in J/mol."""
        return self.phases[0](ca, T), self.phases[1](cb, T)

    # -- tie lines ----------------------------------------------------------------------
    def _seeds(self, temps):
        """Two-phase compositions from pycalphad equilibrium on a composition grid."""
        xs = np.arange(0.0025, 1.0, 0.005)
        pa, pb = self.phase_names
        with self.system._lock:
            eq = equilibrium(self.system.db, self.comps, [pa, pb],
                             {v.X(self.B): xs, v.T: np.atleast_1d(temps),
                              v.P: P_DEFAULT, v.N: 1})
        phase = eq.Phase.values[0, 0]                  # (T, x, vertex)
        X = eq.X.sel(component=self.B).values[0, 0]
        seeds = []
        for i in range(len(np.atleast_1d(temps))):
            best = None
            for j, xj in enumerate(xs):
                names = list(phase[i, j])
                if pa in names and pb in names:
                    d = abs(xj - self.c_alloy) if self.c_alloy is not None else 0.0
                    if best is None or d < best[0]:
                        best = (d, X[i, j, names.index(pa)], X[i, j, names.index(pb)])
                    if self.c_alloy is None:
                        break
            seeds.append(None if best is None else (best[1], best[2]))
        return seeds

    def _common_tangent(self, T, a, b, energies, tol=1e-10, maxit=100):
        fa, fb = energies
        for _ in range(maxit):
            Ga, dA, d2A = (float(q) for q in fa(a, T))
            Gb, dB, d2B = (float(q) for q in fb(b, T))
            Rv = np.array([dA - dB, (Ga - a * dA) - (Gb - b * dB)])
            if abs(Rv[0]) < tol and abs(Rv[1]) < tol:
                break
            J = np.array([[d2A, -d2B], [-a * d2A, b * d2B]])
            step = np.linalg.solve(J, -Rv)
            step = np.clip(step, -0.01, 0.01)
            a = float(np.clip(a + step[0], 1e-9, 1 - 1e-9))
            b = float(np.clip(b + step[1], 1e-9, 1 - 1e-9))
        return a, b

    def tie_lines(self, temps, energies=None):
        """Equilibrium compositions of the two phases vs T, plus what a KKS code
        needs: mu_eq, G at the tie line, curvatures A = d2G/dc2, k = c_a/c_b.
        Uses the reference (exact/tabulated) energies, seeded by pycalphad."""
        temps = [float(t) for t in np.atleast_1d(temps)]
        energies = energies or self.reference
        todo = [t for t in temps if (t, id(energies)) not in self._tie_cache]
        if todo:
            for t, s in zip(todo, self._seeds(todo)):
                row = None
                if s is not None:
                    a, b = self._common_tangent(t, s[0], s[1], energies)
                    Ga, mu, Aa = (float(q) for q in energies[0](a, t))
                    Gb, _, Ab = (float(q) for q in energies[1](b, t))
                    row = {"T": t, "c_a": a, "c_b": b, "k": a / b if b else None,
                           "mu": mu, "G_a": Ga, "G_b": Gb, "A_a": Aa, "A_b": Ab,
                           "c_a_pycalphad": float(s[0]), "c_b_pycalphad": float(s[1]),
                           "consistent": bool(abs(a - s[0]) < 1e-3 and abs(b - s[1]) < 1e-3)}
                self._tie_cache[(t, id(energies))] = row
        return [self._tie_cache[(t, id(energies))] for t in temps]

    def tie_line(self, T):
        return self.tie_lines([T])[0]

    # -- the KKS constraint ----------------------------------------------------------------
    def split(self, c, phi=None, T=None, h=None, ca0=None, cb0=None,
              tol=1e-10, maxit=60, return_info=False):
        """Solve, at every point,
             c = h c_a + (1 - h) c_b    and    dG_a/dc(c_a) = dG_b/dc(c_b) = mu.
        Give phi (h is computed) or h directly. T: scalar or array like c.
        ca0, cb0: initial guesses (pass the previous time step's values for speed).
        Returns c_a, c_b, mu [J/mol] (and an info dict if return_info)."""
        if T is None:
            raise CalphadError("split needs T")
        if h is None:
            if phi is None:
                raise CalphadError("give phi or h")
            h = self.h(np.asarray(phi, dtype=float))
        c, T = _bcast(c, T)
        h = np.broadcast_to(np.asarray(h, dtype=float), c.shape)

        if self.mode == "parabolic":
            pa, pb = self.phases
            ea, Aa, mu_eq, _ = pa.coefficients(T)
            eb, Ab, _, _ = pb.coefficients(T)
            r = Aa / Ab
            ca = (c - (1 - h) * (eb - r * ea)) / (h + (1 - h) * r)
            cb = eb + r * (ca - ea)
            mu = mu_eq + Aa * (ca - ea)
            info = {"iterations": 0, "max_residual": 0.0}
            return (ca, cb, mu, info) if return_info else (ca, cb, mu)

        ca = np.array(c if ca0 is None else np.broadcast_to(ca0, c.shape), dtype=float)
        cb = np.array(c if cb0 is None else np.broadcast_to(cb0, c.shape), dtype=float)
        lo, hi = 1e-10, 1 - 1e-10
        fa, fb = self.phases
        for it in range(1, maxit + 1):
            _, dA, d2A = fa(ca, T)
            _, dB, d2B = fb(cb, T)
            R1 = h * ca + (1 - h) * cb - c
            R2 = dA - dB
            scale = np.maximum(np.abs(d2A), np.abs(d2B))
            err = max(np.max(np.abs(R1)), np.max(np.abs(R2) / scale))
            if err < tol:
                break
            det = -h * d2B - (1 - h) * d2A               # J = [[h, 1-h], [d2A, -d2B]]
            dca = (R1 * d2B + (1 - h) * R2) / det
            dcb = (-h * R2 + d2A * R1) / det
            lim = 0.1                                    # damp large steps
            m = np.maximum(np.maximum(np.abs(dca), np.abs(dcb)) / lim, 1.0)
            ca = np.clip(ca + dca / m, lo, hi)
            cb = np.clip(cb + dcb / m, lo, hi)
        _, mu, _ = fa(ca, T)
        info = {"iterations": it, "max_residual": float(err)}
        if err > 1e3 * tol:
            info["warning"] = "KKS split did not fully converge"
        return (ca, cb, mu, info) if return_info else (ca, cb, mu)

    # -- quantities for the evolution equations --------------------------------------
    def _per_volume(self, val, per_volume):
        if not per_volume:
            return val
        if not self.Vm:
            raise CalphadError("Set Vm (m^3/mol) for per-volume quantities.")
        return val / self.Vm

    def driving_force(self, ca, cb, mu, T, per_volume=True):
        """f_a(c_a) - f_b(c_b) - mu (c_a - c_b). Negative favours phase a (phi -> 1).
        J/m^3 by default (J/mol with per_volume=False). Enters the phase-field
        equation as  - h'(phi) * driving_force."""
        (Ga, _, _), (Gb, _, _) = self.energies(ca, cb, T)
        return self._per_volume(Ga - Gb - mu * (ca - cb), per_volume)

    def free_energy_density(self, c, phi, T, per_volume=True):
        """Chemical part h f_a(c_a) + (1-h) f_b(c_b) (no double well, no gradient)."""
        ca, cb, _ = self.split(c, phi, T)
        h = self.h(np.asarray(phi, dtype=float))
        (Ga, _, _), (Gb, _, _) = self.energies(ca, cb, T)
        return self._per_volume(h * Ga + (1 - h) * Gb, per_volume)

    def dmu_dc(self, ca, cb, phi, T):
        """d(mu)/dc at fixed phi = 1 / ( h/G_a'' + (1-h)/G_b'' )  [J/mol].
        For the mu-form of the diffusion equation: dc/dt = div( D/dmu_dc grad mu )."""
        h = self.h(np.asarray(phi, dtype=float))
        (_, _, d2A), (_, _, d2B) = self.energies(ca, cb, T)
        return 1.0 / (h / d2A + (1 - h) / d2B)

    # -- export for codes in other languages ----------------------------------------------
    def parameters(self, temps, c_grid=None, include_tables=False):
        """Everything a KKS solver needs, as plain data (JSON-friendly)."""
        rows = self.tie_lines(temps)
        out = {
            "solute": self.B, "solvent": self.A,
            "phase_a": self.phase_names[0], "phase_b": self.phase_names[1],
            "mode": self.mode, "Vm": self.Vm,
            "convention": "c = x_solute; phase_a at phi=1; energies J/mol; "
                          "f_parabolic = G_eq + mu_eq (c - c_eq) + A/2 (c - c_eq)^2",
            "tie_lines": rows,
        }
        if include_tables:
            c_grid = np.linspace(1e-3, 1 - 1e-3, 199) if c_grid is None else np.asarray(c_grid)
            tabs = {}
            for ph in self.phases:
                G, d1, d2 = [], [], []
                for t in np.atleast_1d(temps):
                    g0, g1, g2 = ph(c_grid, float(t))
                    G.append(g0.tolist()); d1.append(g1.tolist()); d2.append(g2.tolist())
                tabs[ph.name] = {"G": G, "dG": d1, "d2G": d2}
            out["tables"] = {"T": [float(t) for t in np.atleast_1d(temps)],
                             "c": c_grid.tolist(), "phases": tabs}
        return out

    def export_parabolic_csv(self, path, temps):
        """CSV of tie lines and parabolic coefficients vs T (read by cpp/kks_parabolic.hpp)."""
        rows = [r for r in self.tie_lines(temps) if r is not None]
        with open(path, "w") as fh:
            fh.write(f"# solute={self.B} solvent={self.A} phase_a={self.phase_names[0]} "
                     f"phase_b={self.phase_names[1]} Vm={self.Vm}\n")
            fh.write("T_K,c_a_eq,c_b_eq,mu_eq,G_a_eq,G_b_eq,A_a,A_b\n")
            for r in rows:
                fh.write(f"{r['T']:.6f},{r['c_a']:.12g},{r['c_b']:.12g},{r['mu']:.12g},"
                         f"{r['G_a']:.12g},{r['G_b']:.12g},{r['A_a']:.12g},{r['A_b']:.12g}\n")
        return len(rows)

    def export_json(self, path, temps, include_tables=True, c_grid=None):
        with open(path, "w") as fh:
            json.dump(self.parameters(temps, c_grid, include_tables), fh)
