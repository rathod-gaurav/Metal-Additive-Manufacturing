"""
calphad_client.py
=================
Minimal Python client for calphad_server.py. Standard library only, so you can
drop this one file into any project.

    from calphad_client import CalphadClient
    c = CalphadClient("http://127.0.0.1:8000")
    r = c.equilibrium({"NI": 0.33}, T=1500, balance="CU")
    for ph in r["results"]["phases"]:
        print(ph["phase"], ph["fraction"], ph["composition"])
"""

import json
import urllib.error
import urllib.request


class CalphadServerError(RuntimeError):
    pass


class CalphadClient:
    def __init__(self, url="http://127.0.0.1:8000", db=None, timeout=600):
        self.url = url.rstrip("/")
        self.db = db          # default database name for every call
        self.timeout = timeout

    # -- transport ------------------------------------------------------------
    def _request(self, method, path, body=None):
        data = None if body is None else json.dumps(body).encode()
        req = urllib.request.Request(self.url + path, data=data, method=method,
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as exc:
            try:
                detail = json.loads(exc.read()).get("detail")
            except Exception:
                detail = exc.reason
            raise CalphadServerError(f"HTTP {exc.code}: {detail}") from None

    def _body(self, db, **kw):
        body = {k: val for k, val in kw.items() if val is not None}
        if db or self.db:
            body["db"] = db or self.db
        return body

    # -- databases -------------------------------------------------------------
    def health(self):
        return self._request("GET", "/health")

    def databases(self):
        return self._request("GET", "/databases")

    def info(self, db=None):
        return self._request("GET", f"/databases/{db or self.db}")

    def load_database(self, path=None, tdb_text=None, name=None):
        """Load a TDB on the server, by server-side path or by file contents."""
        return self._request("POST", "/databases",
                             {"path": path, "tdb_text": tdb_text, "name": name})

    def upload_tdb(self, local_path, name=None):
        """Send a local TDB file's contents to the server."""
        with open(local_path, encoding="utf-8", errors="replace") as fh:
            return self.load_database(tdb_text=fh.read(), name=name)

    # -- the three APIs ----------------------------------------------------------
    def equilibrium(self, composition, T, P=None, phases=None, units=None,
                    balance=None, verify=None, db=None):
        """API 1. T may be a number or a list."""
        return self._request("POST", "/equilibrium", self._body(
            db, composition=composition, T=T, P=P, phases=phases,
            units=units, balance=balance, verify=verify))

    def phase_energies(self, composition, T, P=None, phases=None, units=None,
                       balance=None, db=None):
        """API 2."""
        return self._request("POST", "/phase_energies", self._body(
            db, composition=composition, T=T, P=P, phases=phases,
            units=units, balance=balance))

    def energy_surface(self, T, P=None, elements=None, phases=None, mode=None,
                       step=None, pdens=None, max_points=None, db=None):
        """API 3."""
        return self._request("POST", "/energy_surface", self._body(
            db, T=T, P=P, elements=elements, phases=phases, mode=mode,
            step=step, pdens=pdens, max_points=max_points))

    # -- KKS ------------------------------------------------------------------------
    def kks_parameters(self, solute, phase_a, phase_b, T, mode=None, c_alloy=None,
                       Vm=None, solvent=None, branches=None, include_tables=None,
                       n_c=None, db=None):
        """Tie lines, parabolic coefficients and (optionally) G tables for a
        binary two-phase KKS model. phase_a is the phase at phi = 1."""
        return self._request("POST", "/kks/parameters", self._body(
            db, solute=solute, phase_a=phase_a, phase_b=phase_b, T=list(T), mode=mode,
            c_alloy=c_alloy, Vm=Vm, solvent=solvent, branches=branches,
            include_tables=include_tables, n_c=n_c))


if __name__ == "__main__":
    import sys
    c = CalphadClient(sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8000")
    print(json.dumps(c.health(), indent=2))
