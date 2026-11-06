"""
calphad_core.py
===============
A thin, reusable layer over pycalphad. Load any TDB file and ask three
questions of it:

  1. equilibrium()     fixed alloy composition + T  ->  stable phases,
                       phase fractions, phase compositions, site fractions,
                       chemical potentials.
  2. phase_energies()  fixed alloy composition + T  ->  molar Gibbs energy of
                       EACH phase at that composition (whether stable or not),
                       plus that phase's chemical potentials there.
  3. energy_surface()  T  ->  G(x) for every phase over the whole composition
                       space (curves for binaries, surfaces for ternaries).

This module has no web code in it. Import it directly from Python, or run
calphad_server.py to expose the same functions over HTTP/JSON so any
language can call them.

Conventions
-----------
* Compositions are mole fractions of elements unless `units` says otherwise.
* All energies are J per mole of atoms, referenced to the database's
  reference state (usually SER: pure elements, stable structure, 298.15 K).
* Phase fractions are fractions of the total moles of atoms.
* P is in Pa, T in K.
"""

import itertools
import math
import os
import threading
import warnings

import numpy as np
from pycalphad import Database, Model, Workspace, calculate, equilibrium, variables as v
from pycalphad.property_framework.metaproperties import DormantPhase, IsolatedPhase
from pycalphad.core.utils import filter_phases

warnings.filterwarnings("ignore", category=UserWarning)

P_DEFAULT = 101325.0
UNITS = ("mole_fraction", "atomic_percent", "weight_fraction", "weight_percent")


class CalphadError(ValueError):
    """Bad input (unknown element/phase, composition not summing to 1, ...)."""


def _num(x):
    """numpy scalar -> float, NaN/inf -> None (JSON has no NaN)."""
    x = float(x)
    return x if math.isfinite(x) else None


class CalphadSystem:
    """One loaded thermodynamic database."""

    def __init__(self, db, name="db"):
        self.db = db
        self.name = name
        self.elements = sorted(e for e in db.elements if e not in ("VA", "/-"))
        self.phases = sorted(db.phases.keys())
        self._models = {}
        self._lock = threading.RLock()  # pycalphad calls are serialised per database

    # -- construction ------------------------------------------------------
    @classmethod
    def from_file(cls, path, name=None):
        if not os.path.isfile(path):
            raise CalphadError(f"TDB file not found: {path}")
        name = name or os.path.splitext(os.path.basename(path))[0]
        return cls(Database(path), name)

    @classmethod
    def from_text(cls, text, name="db"):
        return cls(Database(text), name)

    # -- database description ----------------------------------------------
    def info(self):
        phases = {}
        for ph in self.phases:
            p = self.db.phases[ph]
            hints = p.model_hints or {}
            entry = {
                "site_ratios": [float(r) for r in p.sublattices],
                "constituents": [sorted(s.name for s in sl) for sl in p.constituents],
            }
            if "ordered_phase" in hints and "disordered_phase" in hints:
                entry["order_disorder_pair"] = {
                    "ordered": hints["ordered_phase"],
                    "disordered": hints["disordered_phase"]}
            phases[ph] = entry
        masses = {}
        for el in self.elements:
            m = self.db.refstates.get(el, {}).get("mass")
            masses[el] = _num(m) if m is not None else None
        return {"name": self.name, "elements": self.elements,
                "atomic_masses": masses, "phases": phases}

    # -- helpers --------------------------------------------------------------
    def _components(self, elements):
        comps = sorted(elements)
        if "VA" in self.db.elements:
            comps.append("VA")
        return comps

    def _check_elements(self, elements):
        bad = [e for e in elements if e not in self.elements]
        if bad:
            raise CalphadError(f"Unknown element(s) {bad}. Database has {self.elements}.")

    def _parse_composition(self, composition, units="mole_fraction", balance=None):
        """Return {element: mole fraction} that sums to 1.

        composition : dict like {"NI": 0.3} or {"CU": 70, "NI": 30}
        units       : mole_fraction | atomic_percent | weight_fraction | weight_percent
        balance     : element that takes up the remainder (optional). If omitted,
                      the given values must already sum to 1 (or 100).
        """
        if units not in UNITS:
            raise CalphadError(f"units must be one of {UNITS}")
        comp = {str(k).upper(): float(val) for k, val in composition.items()}
        if balance:
            balance = balance.upper()
            comp.pop(balance, None)
        self._check_elements(list(comp) + ([balance] if balance else []))
        if any(val < 0 for val in comp.values()):
            raise CalphadError("Composition values must be non-negative.")
        scale = 100.0 if units.endswith("percent") else 1.0
        comp = {k: val / scale for k, val in comp.items()}
        if balance:
            rest = 1.0 - sum(comp.values())
            if rest < -1e-9:
                raise CalphadError("Solute amounts exceed 100%.")
            comp[balance] = max(rest, 0.0)
        total = sum(comp.values())
        if abs(total - 1.0) > 1e-6:
            raise CalphadError(f"Composition sums to {total*scale:g}, not {scale:g}. "
                               "Give all elements or name a 'balance' element.")
        if units.startswith("weight"):
            masses = self.info()["atomic_masses"]
            if any(masses[e] is None for e in comp):
                raise CalphadError("Database lacks atomic masses; use mole fractions.")
            moles = {e: w / masses[e] for e, w in comp.items()}
            tot = sum(moles.values())
            comp = {e: n / tot for e, n in moles.items()}
        comp = {e: x for e, x in comp.items()}
        if len([e for e, x in comp.items() if x > 0]) < 1:
            raise CalphadError("Empty composition.")
        return comp

    def _conditions(self, x, T, P):
        elems = sorted(x)
        dependent = max(elems, key=lambda e: x[e])  # any one element is implied
        conds = {v.X(e): x[e] for e in elems if e != dependent}
        conds.update({v.T: T, v.P: P, v.N: 1})
        return conds

    def _phase_list(self, comps, phases):
        if phases:
            phases = [p.upper() for p in phases]
            bad = [p for p in phases if p not in self.db.phases]
            if bad:
                raise CalphadError(f"Unknown phase(s) {bad}. Database has {self.phases}.")
            return filter_phases(self.db, comps, phases)
        # Default: every phase that can form from these elements. For
        # order/disorder pairs pycalphad keeps the ordered phase, which
        # describes both states (see 'state' in the results).
        return filter_phases(self.db, comps)

    def _model(self, comps, phase):
        key = (tuple(comps), phase)
        if key not in self._models:
            self._models[key] = Model(self.db, comps, phase)
        return self._models[key]

    def _models_for(self, comps, phases):
        return {ph: self._model(comps, ph) for ph in phases}

    def _site_fraction_dict(self, model, y):
        out = [dict() for _ in model.site_ratios]
        for sf, yy in zip(model.site_fractions, y):
            out[sf.sublattice_index][sf.species.name] = _num(yy)
        return out

    def _state(self, phase, sublattices):
        """Label an order/disorder phase as ordered or disordered from its
        site fractions (identical sublattices = disordered)."""
        hints = self.db.phases[phase].model_hints or {}
        # compare only the substitutional sublattices (skip vacancy-only ones)
        subs = [sl for sl in sublattices if set(sl) - {"VA"}]
        if hints.get("ordered_phase") != phase or len(subs) < 2:
            return None, None
        els = set().union(*subs)
        spread = max(abs((subs[i].get(e) or 0) - (subs[0].get(e) or 0))
                     for i in range(1, len(subs)) for e in els)
        if spread < 1e-3:
            return "disordered", hints.get("disordered_phase")
        return "ordered", phase

    def _mass_fractions(self, xdict):
        masses = self.info()["atomic_masses"]
        if any(masses.get(e) is None for e in xdict):
            return None
        m = {e: (x or 0) * masses[e] for e, x in xdict.items()}
        tot = sum(m.values())
        return {e: (val / tot if tot > 0 else None) for e, val in m.items()}

    # =======================================================================
    # One phase on its own: robust homogeneous minimum and driving force
    # =======================================================================
    # pycalphad's local minimisers can stop at the metastable state of an
    # order/disorder phase (seen: Ni-25Al gamma-prime at 1600 K returned
    # disordered, 355 J/mol too high). So we also refine from several explicit
    # starting configurations (disordered, and each element enriched on each
    # sublattice) and keep the lowest.
    def _starting_points(self, model, x):
        subl = [sorted(sp.name for sp in sl) for sl in model.constituents]
        ratios = [float(r) for r in model.site_ratios]

        def base():
            y = {}
            for s, names in enumerate(subl):
                real = [n for n in names if n != "VA"]
                if not real:
                    y[(s, "VA")] = 1.0
                    continue
                tot = sum(x.get(n, 0.0) for n in real)
                for n in real:
                    y[(s, n)] = x.get(n, 0.0) / tot if tot > 0 else 1.0 / len(real)
                if "VA" in names:
                    y[(s, "VA")] = 1e-4
            return y

        starts = [base()]
        mixing = [s for s, names in enumerate(subl) if len([n for n in names if n != "VA"]) >= 2]
        if len(mixing) >= 2:
            for s in mixing:
                real = [n for n in subl[s] if n != "VA"]
                for el in real:
                    if x.get(el, 0.0) <= 0:
                        continue
                    y = base()
                    for n in real:
                        y[(s, n)] = 0.9 if n == el else 0.1 / (len(real) - 1)
                    a_s = ratios[s]
                    a_o = sum(ratios[t] for t in mixing if t != s)
                    need = min(max((x[el] * (a_s + a_o) - 0.9 * a_s) / a_o, 1e-3), 0.999)
                    for t in mixing:
                        if t == s:
                            continue
                        rt = [n for n in subl[t] if n != "VA"]
                        if el not in rt:
                            continue
                        rest = [n for n in rt if n != el]
                        tot = sum(y[(t, n)] for n in rest) or 1.0
                        y[(t, el)] = need
                        for n in rest:
                            y[(t, n)] = (1 - need) * y[(t, n)] / tot
                    starts.append(y)
        out = []
        for y in starts:
            raw = [max(y.get((sf.sublattice_index, sf.species.name), 1e-12), 1e-12)
                   for sf in model.site_fractions]
            sums = {}
            for sf, val in zip(model.site_fractions, raw):
                sums[sf.sublattice_index] = sums.get(sf.sublattice_index, 0.0) + val
            out.append(np.array([val / sums[sf.sublattice_index]
                                 for sf, val in zip(model.site_fractions, raw)]))
        return out

    def _refine(self, comps, phase, x, temps, P, y0):
        """Minimise one phase as a single composition set at composition x,
        starting from site fractions y0. Vectorised over temps.
        Returns GM (T,), MU (T, elements), Y (T, n_sitefractions); NaN if failed."""
        model = self._model(comps, phase)
        temps = np.atleast_1d(np.asarray(temps, dtype=float))
        conds = self._conditions(x, list(temps), P)
        eq = equilibrium(self.db, comps, [phase], conds, model={phase: model},
                         calc_opts={"points": {phase: np.atleast_2d(y0)}})
        e = eq.isel({d: 0 for d in eq.sizes if d in ("N", "P") or d.startswith("X_")})
        nsf = len(model.site_fractions)
        GM = np.full(len(temps), np.nan)
        MU = np.full((len(temps), e.sizes["component"]), np.nan)
        Y = np.full((len(temps), nsf), np.nan)
        for i in range(len(temps)):
            names = [str(p) for p in e.Phase.values[i]]
            sets = [k for k, n in enumerate(names) if n == phase]
            if len(sets) == 1 and np.isfinite(e.GM.values[i]):
                GM[i] = e.GM.values[i]
                MU[i] = e.MU.values[i]
                Y[i] = e.Y.values[i, sets[0], :nsf]
        return GM, MU, Y, [str(c) for c in e.component.values]

    def _homogeneous_minimum(self, comps, phase, x, T, P):
        """Lowest single-composition-set state of `phase` at composition x."""
        model = self._model(comps, phase)
        best = None
        for y0 in self._starting_points(model, x):
            try:
                GM, MU, Y, els = self._refine(comps, phase, x, [T], P, y0)
            except Exception:
                continue
            if np.isfinite(GM[0]) and (best is None or GM[0] < best[0] - 1e-9):
                best = (float(GM[0]), dict(zip(els, MU[0])), Y[0])
        return best

    def _driving_force(self, wks, comps, phase, mu, T, P):
        """max(pycalphad DormantPhase, sampled estimate) in J/mol.
        DormantPhase is a local search and was seen to come out ~11 J/mol low
        for gamma-prime; the sampled estimate max(sum x_i mu_i - G) guards it."""
        vals = []
        try:
            vals.append(_num(np.ravel(wks.get(DormantPhase(phase, wks=wks).driving_force))[0]))
        except Exception:
            pass
        try:
            model = self._model(comps, phase)
            res = calculate(self.db, comps, phase, T=T, P=P, N=1, pdens=2000,
                            output="GM", model={phase: model})
            els = [str(c) for c in res.component.values]
            X = res.X.values.reshape(-1, len(els))
            G = res.GM.values.reshape(-1)
            m = np.array([mu[e] for e in els])
            good = np.isfinite(G)
            if good.any():
                vals.append(_num(np.max(X[good] @ m - G[good])))
        except Exception:
            pass
        vals = [v_ for v_ in vals if v_ is not None]
        return max(vals) if vals else None

    # =======================================================================
    # Equilibrium solve with a global-minimum check
    # =======================================================================
    # pycalphad starts its minimiser from a sampled grid of configurations. If
    # that grid misses a phase's low-energy region, the solver can settle in a
    # local minimum (we saw LIQUID+L12 reported where L12+B2 is lower). The check:
    # every absent phase must have driving force <= DF_TOL against the solution's
    # chemical potentials. If not, re-solve with denser sampling.
    DF_TOL = 1.0                       # J/mol
    PDENS_LADDER = (None, 2000, 10000)  # None = pycalphad default

    def _solve(self, x, comps, phase_list, T, P, verify=True):
        conds = self._conditions(x, T, P)
        models = self._models_for(comps, phase_list)
        ladder = self.PDENS_LADDER if verify else (None,)
        for pd in ladder:
            opts = {} if pd is None else {"pdens": pd}
            wks = Workspace(self.db, comps, phase_list, conds, models=models, calc_opts=opts)
            present = {str(p) for p in np.ravel(wks.eq.Phase) if str(p)}
            dfs = {}
            if verify:
                els = [str(c) for c in wks.eq.coords["component"]]
                mu = {el: float(m) for el, m in zip(els, np.ravel(wks.eq.MU))}
                for ph in phase_list:
                    if ph in present:
                        continue
                    dfs[ph] = self._driving_force(wks, comps, ph, mu, T, P)
            worst = max([d for d in dfs.values() if d is not None], default=None)
            check = {"verified": bool(verify and (worst is None or worst <= self.DF_TOL)),
                     "sampling_density": pd or "default",
                     "driving_forces_absent_phases": dfs if verify else None}
            if not verify or check["verified"]:
                break
        return wks, models, check

    def _parse(self, e, models):
        """One equilibrium point (pycalphad LightDataset) -> dict."""
        phases = [str(p) for p in np.ravel(e.Phase)]
        nv = len(phases)
        NP = np.ravel(e.NP)
        X = np.asarray(e.X).reshape(nv, -1)
        Y = np.asarray(e.Y).reshape(nv, -1)
        el_order = [str(c) for c in e.coords["component"]]
        mu = {el: _num(m) for el, m in zip(el_order, np.ravel(e.MU))}
        stable, seen = [], {}
        for k, ph in enumerate(phases):
            if not ph:
                continue
            seen[ph] = seen.get(ph, 0) + 1
            xph = {el: _num(val) for el, val in zip(el_order, X[k])}
            model = models[ph]
            sl = self._site_fraction_dict(model, Y[k][:len(model.site_fractions)])
            state, shown_as = self._state(ph, sl)
            # every stable phase touches the equilibrium tangent plane
            gm = sum((xph[el] or 0) * (mu[el] or 0) for el in el_order)
            stable.append({
                "phase": ph,
                "composition_set": seen[ph],
                "state": state,
                "equivalent_to": shown_as,
                "fraction": _num(NP[k]),
                "composition": xph,
                "composition_wt": self._mass_fractions(xph),
                "GM": _num(gm),
                "site_fractions": sl,
            })
        stable.sort(key=lambda d: -(d["fraction"] or 0))
        return {"converged": bool(stable) and all(d["fraction"] is not None for d in stable),
                "GM_system": _num(np.ravel(e.GM)[0]),
                "chemical_potentials": mu,
                "phases": stable}

    # =======================================================================
    # API 1: equilibrium
    # =======================================================================
    def equilibrium(self, composition, T, P=P_DEFAULT, phases=None,
                    units="mole_fraction", balance=None, verify=True):
        """Stable phases at fixed overall composition.

        T may be a number or a list (one result per T).
        verify=True checks that no absent phase has a positive driving force
        and re-solves with denser sampling if one does (see _solve).
        """
        x = self._parse_composition(composition, units, balance)
        comps = self._components(x)
        phase_list = self._phase_list(comps, phases)
        temps = [float(t) for t in np.atleast_1d(T)]
        results = []
        with self._lock:
            for t in temps:
                wks, models, check = self._solve(x, comps, phase_list, t, float(P), verify)
                res = {"T": t, "P": float(P)}
                res.update(self._parse(wks.eq, models))
                res["stability_check"] = check
                results.append(res)
        return {"alloy": x, "units_in": units, "candidate_phases": phase_list,
                "results": results if isinstance(T, (list, tuple, np.ndarray)) else results[0]}

    # =======================================================================
    # API 2: Gibbs energy of every phase at the alloy composition
    # =======================================================================
    def phase_energies(self, composition, T, P=P_DEFAULT, phases=None,
                       units="mole_fraction", balance=None):
        """G of each phase at the overall alloy composition, stable or not.

        For each phase, pycalphad's IsolatedPhase minimises that phase ALONE
        as a single homogeneous composition set at the alloy composition
        (sublattice occupancies / ordering are relaxed, composition is fixed).
        Returned per phase:
          GM                     J/mol of atoms
          chemical_potentials    that phase's tangent plane at this composition.
                                 Differences MU(i) - MU(dependent) are the
                                 diffusion potentials a KKS code uses.
          GM_above_equilibrium   GM - G of the equilibrium phase mixture (>= 0)
          driving_force          J/mol, against the equilibrium chemical
                                 potentials: 0 for stable phases, < 0 for phases
                                 that would not form (pycalphad DormantPhase).
          splits_if_allowed      present when the phase on its own would rather
                                 separate into two compositions (e.g. ordered +
                                 disordered, or a miscibility gap).
        """
        x = self._parse_composition(composition, units, balance)
        comps = self._components(x)
        elements = [c for c in comps if c != "VA"]
        phase_list = self._phase_list(comps, phases)
        T, P = float(T), float(P)
        conds = self._conditions(x, T, P)

        out = []
        with self._lock:
            wks, models, check = self._solve(x, comps, phase_list, T, P)
            ref = self._parse(wks.eq, models)
            g_eq = ref["GM_system"]
            pd = check["sampling_density"]
            opts = {} if pd == "default" else {"pdens": pd}
            for ph in phase_list:
                model = models[ph]
                entry = {"phase": ph}
                try:
                    ip = IsolatedPhase(ph, wks=wks)
                    props = (["GM"] + [f"MU({el})" for el in elements]
                             + [str(sf) for sf in model.site_fractions])
                    vals = [_num(np.ravel(val)[0]) for val in wks.get(*[ip(q) for q in props])]
                except Exception as exc:  # pycalphad could not converge this phase
                    vals = [None]
                    entry["note"] = f"isolated-phase calculation failed: {exc}"
                # second opinion: refine from explicit ordered/disordered starts
                robust = self._homogeneous_minimum(comps, ph, x, T, P)
                if robust is not None and (vals[0] is None or robust[0] < vals[0] - 1e-6):
                    entry.pop("note", None)
                    vals = ([robust[0]] + [_num(robust[1].get(el)) for el in elements]
                            + [_num(yy) for yy in robust[2]])
                gm = vals[0]
                entry["GM"] = gm
                if gm is None:
                    entry.setdefault("note", "phase cannot attain this composition "
                                             "(e.g. a stoichiometric compound)")
                    out.append(entry)
                    continue
                n = len(elements)
                entry["chemical_potentials"] = dict(zip(elements, vals[1:1 + n]))
                sl = self._site_fraction_dict(model, vals[1 + n:])
                state, shown_as = self._state(ph, sl)
                entry.update(state=state, equivalent_to=shown_as, site_fractions=sl)
                entry["GM_above_equilibrium"] = _num(gm - g_eq) if g_eq is not None else None
                if ph in {p["phase"] for p in ref["phases"]}:
                    entry["driving_force"] = 0.0          # stable: on the tangent plane
                else:
                    entry["driving_force"] = self._driving_force(
                        wks, comps, ph, ref["chemical_potentials"], T, P)

                # Would this phase, on its own, rather split into two compositions?
                eq1 = equilibrium(self.db, comps, [ph], conds, model={ph: model},
                                  calc_opts=opts)
                e = eq1.isel({d: 0 for d in eq1.sizes
                              if d not in ("vertex", "component", "internal_dof")})
                sets = [k for k in range(e.sizes["vertex"]) if str(e.Phase.values[k]) == ph]
                g_split = _num(e.GM.values)
                if len(sets) > 1 and g_split is not None and g_split < gm - 1e-6:
                    el_order = [str(c) for c in e.component.values]
                    entry["splits_if_allowed"] = {
                        "GM": g_split,
                        "composition_sets": [{
                            "fraction": _num(e.NP.values[k]),
                            "composition": {el: _num(val) for el, val in
                                            zip(el_order, e.X.values[k])},
                        } for k in sets]}
                out.append(entry)
        out.sort(key=lambda d: (d["GM"] is None, d["GM"] if d["GM"] is not None else 0))
        return {"alloy": x, "T": T, "P": P, "GM_equilibrium": g_eq,
                "stable_phases": [p["phase"] + (f"#{p['composition_set']}"
                                                if p["composition_set"] > 1 else "")
                                  for p in ref["phases"]],
                "stability_check": check,
                "phases": out}

    # =======================================================================
    # API 3: Gibbs energy surfaces
    # =======================================================================
    def _mixing_sublattice(self, model, elements):
        """Index of the single sublattice that holds all the elements, if the
        phase is a plain substitutional solution (other sublattices VA only)."""
        idx = None
        for s, sl in enumerate(model.constituents):
            names = {sp.name for sp in sl}
            if names - {"VA"}:
                if idx is not None or "VA" in names:
                    return None
                idx = s
        return idx

    def energy_surface(self, T, P=P_DEFAULT, elements=None, phases=None,
                       step=0.02, mode="grid", pdens=500, max_points=200000):
        """G(x) of every phase at temperature T.

        mode="grid"   : values on a regular composition grid with spacing `step`.
                        Exact for substitutional solutions; for phases with
                        sublattices it is the lowest sampled G in each grid cell
                        (flagged "approximate": true; raise pdens to refine).
        mode="points" : the raw sampled points from pycalphad.calculate,
                        including every internal configuration.
        """
        elements = sorted(e.upper() for e in (elements or self.elements))
        self._check_elements(elements)
        comps = self._components(elements)
        phase_list = self._phase_list(comps, phases)
        n = len(elements)
        T, P = float(T), float(P)

        grid = None
        if mode == "grid":
            K = int(round(1.0 / step))
            count = math.comb(K + n - 1, n - 1)
            if count > max_points:
                raise CalphadError(f"{count} grid points for {n} elements at step={step}; "
                                   f"increase step or max_points ({max_points}).")
            grid = np.array([c for c in itertools.product(range(K + 1), repeat=n - 1)
                             if sum(c) <= K], dtype=float)
            grid = np.column_stack([grid, K - grid.sum(axis=1)]) / K  # rows sum to 1
        elif mode != "points":
            raise CalphadError("mode must be 'grid' or 'points'")

        surfaces = []
        for ph in phase_list:
            with self._lock:
                model = self._model(comps, ph)
                mix = self._mixing_sublattice(model, elements) if mode == "grid" else None
                if mode == "grid" and mix is not None:
                    # exact: site fractions on the mixing sublattice = mole fractions
                    present = {sp.name for sp in model.constituents[mix]}
                    missing = [i for i, e in enumerate(elements) if e not in present]
                    ok = (np.all(grid[:, missing] < 1e-12, axis=1) if missing
                          else np.ones(len(grid), bool))
                    pts_x = grid[ok]
                    cols = []
                    for sf in model.site_fractions:
                        if sf.sublattice_index == mix:
                            cols.append(pts_x[:, elements.index(sf.species.name)])
                        else:
                            cols.append(np.ones(len(pts_x)))  # VA-only sublattice
                    pts = np.clip(np.column_stack(cols), 1e-14, 1.0)
                    res = calculate(self.db, comps, ph, T=T, P=P, N=1, points=pts,
                                    output="GM", model={ph: model})
                    X = res.X.values.reshape(-1, n)
                    G = res.GM.values.reshape(-1)
                    approx = False
                else:
                    res = calculate(self.db, comps, ph, T=T, P=P, N=1, pdens=pdens,
                                    output="GM", model={ph: model})
                    X = res.X.values.reshape(-1, n)
                    G = res.GM.values.reshape(-1)
                    good = np.isfinite(G)
                    X, G = X[good], G[good]
                    # internal degrees of freedom beyond composition -> sampled, not exact
                    internal = len(model.site_fractions) - len(model.site_ratios)
                    approx = mode == "grid" and internal > n - 1
                    if mode == "grid":
                        # keep the lowest-G sample in each grid cell: the phase's
                        # lower envelope over its internal configurations
                        K = int(round(1.0 / step))
                        order = np.argsort(G, kind="stable")
                        X, G = X[order], G[order]
                        keys = np.round(X[:, :-1] * K).astype(int)
                        _, first = np.unique(keys, axis=0, return_index=True)
                        X, G = X[first], G[first]
                    elif len(G) > max_points:
                        keep = np.random.default_rng(0).choice(len(G), max_points, replace=False)
                        X, G = X[keep], G[keep]
            surfaces.append({
                "phase": ph, "n_points": int(len(G)),
                "approximate": bool(approx),
                "X": np.round(X, 10).tolist(),
                "GM": [_num(g) for g in G],
            })
        return {"T": T, "P": P, "elements": elements, "mode": mode,
                "step": step if mode == "grid" else None, "phases": surfaces}
