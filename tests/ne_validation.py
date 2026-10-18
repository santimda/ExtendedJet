#!/usr/bin/env python3
"""Validation of the electron-distribution solver against analytic and semi-analytic solutions.

The solver (edistribution::calcNE_paramo) evolves, for each comoving Lorentz factor, the number flux through the jet
cross-section, F(g, z) = n'(g) * pi r^2 * betaGamma * c, along z:

    dF/dz = - d(gdot_z F)/dg + Q_lin(g, z),     gdot_z = gdot / (betaGamma c) < 0,
    Q_lin = Q_E(E) * m c^2 * pi r^2 * exp(-z/zmax_inj),

with the implicit upwind kernel FP_FinDif_cool. Two groups of tests:

Kernel (FP_FinDif_cool called through the f2py shim tests/ne_kernel_shim.f90, idealised loss and injection terms):
  K1  impulsive injection, pure synchrotron cooling  gdot = -b g^2   (closed form along characteristics)
  K2  impulsive injection, pure adiabatic cooling    gdot = -a g     (closed form)
  K3  continuous injection, synchrotron + adiabatic, steady state    N = int_g^inf Q dg' / |gdot|
  K4  particle-number conservation of K1 (no particle leaves the grid)
  each of K1-K3 is checked at numbins = 1100 and for first-order convergence (550, 1100, 2200 bins)

Full model (run.x on a namelist; N_e only; at the namelist's numbins and at twice that):
  M1  N_e(E, z) from Ne_cuts.dat vs the solution by characteristics with the code's own B(z), Gamma(z), r(z),
      losses (synchrotron + adiabatic), injection normalisation, cutoff E_cut(z) and taper:
          F(g, z) = int_{z0}^{z} Q_lin(g'(z'), z') (dg'/dg) dz',   g'(z') = energy at z' of a particle with g at z.
      This is exact for the transport equation the code solves (no escape, no IC: TLosses = sync + adiabatic).
      Bounded at 5 per cent where injection is active (z <= zmax_inj); first-order convergence required at all z.
  M2  the one-zone (local fast-cooling) limit  N = int_E^inf Q dE' / |Edot|, which M1 must approach where the
      synchrotron cooling length is much shorter than z (reported where t_sync < t_adv/30; shown in the plot).

usage:  python tests/ne_validation.py [--nml configs/m87_reference.nml] [--plots tests/output/ne_validation]
Exit status 1 if any check exceeds its tolerance. Needs run.x and the f2py module (make -C Code; make -C Code f2py).
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
CODE = HERE.parent / "Code"
sys.path.insert(0, str(CODE))
import JetCode as J  # noqa: E402

C_LIGHT = float(J.constants.c)
MEC2 = float(J.constants.mec2)
RESULTS = []


def check(name, value, tol):
    ok = bool(np.isfinite(value) and value <= tol)
    RESULTS.append((name, value, tol, ok))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {value:.3e} (tol {tol:.1e})")
    return ok


# ----------------------------------------------------------------------------------------------------------------
# Kernel tests
# ----------------------------------------------------------------------------------------------------------------
KBUILD = HERE / "output" / "ne_kernel_build"
SHIM = HERE / "ne_kernel_shim.f90"
_KERNEL = None


def kernel():
    """f2py module `ne_kernel` (tests/ne_kernel_shim.f90 + the production sources, the Makefile's FFLAGS).

    FP_FinDif_cool cannot be called through JetCode: f2py cannot size its assumed-shape intent(out) array.
    Rebuilt when the shim or any Code/*.f90 is newer than the module."""
    global _KERNEL
    if _KERNEL is None:
        db = subprocess.run(["make", "-pn", "run.x"], cwd=CODE, capture_output=True, text=True).stdout
        var = {k: v for k, v in re.findall(r"^(FFLAGS|F2PY_SOURCES) = (.*)$", db, flags=re.M)}
        srcs = [CODE / s for s in var["F2PY_SOURCES"].split()]
        KBUILD.mkdir(parents=True, exist_ok=True)
        built = list(KBUILD.glob("ne_kernel*.so"))
        newest = max(p.stat().st_mtime for p in srcs + [SHIM])
        if not built or built[0].stat().st_mtime < newest:
            print("  building the ne_kernel f2py shim (once) ...")
            subprocess.run(["f2py", "-c", *map(str, srcs), str(SHIM), "-m", "ne_kernel", "--backend", "meson",
                            f"--f90flags={var['FFLAGS']} -fopenmp", "-lgomp", "only:", "fp_steps", ":"],
                           cwd=KBUILD, check=True, capture_output=True)
        sys.path.insert(0, str(KBUILD))
        import ne_kernel
        _KERNEL = ne_kernel
    return _KERNEL


def evolve(g, n0, gdot, Q, dt, nsteps):
    return kernel().fp_steps(nsteps, dt, g, n0, gdot, Q, 1e200)


def number(g, n):
    """Total number with the kernel's own cell widths dx = sqrt(dx+ dx-) (the quantity FP_FinDif_cool conserves)."""
    dxp = np.append(np.diff(g), g[-1] - g[-2])
    dxm = np.insert(np.diff(g), 0, g[1] - g[0])
    return np.sum(n * np.sqrt(dxp * dxm))


def powerlaw(g, p, g1, g2):
    return np.where((g >= g1) & (g <= g2), g ** (-p), 0.0)


def k1_sync(nb, nsteps=2000, p=2.2, g1=1e3, g2=1e8, gb=1e5):
    """Impulsive power law cooled by gdot = -b g^2 for t = 1/(b gb): cooling break at gb."""
    g = np.geomspace(2.0, 1e10, nb)
    b, t = 1.0, 1.0 / gb
    n0 = powerlaw(g, p, g1, g2)
    n = evolve(g, n0, -b * g ** 2, np.zeros(nb), t / nsteps, nsteps)
    x = 1.0 - b * g * t
    g0 = np.where(x > 0, g / np.where(x > 0, x, 1), np.inf)
    exact = np.where(x > 0, powerlaw(g0, p, g1, g2) / np.where(x > 0, x, 1) ** 2, 0.0)
    return g, n0, n, exact


def k2_adiabatic(nb=1100, nsteps=2000, p=2.2, g1=1e3, g2=1e8, a=1.0, t=3.0):
    """Impulsive power law cooled by gdot = -a g for t = 3/a (energies shift by e^-3)."""
    g = np.geomspace(2.0, 1e10, nb)
    n0 = powerlaw(g, p, g1, g2)
    n = evolve(g, n0, -a * g, np.zeros(nb), t / nsteps, nsteps)
    exact = powerlaw(g * np.exp(a * t), p, g1, g2) * np.exp(a * t)
    return g, n0, n, exact


def k3_steady(nb=1100, p=2.2, g1=10.0, g2=1e8, a=1.0, b=1e-5):
    """Constant injection Q = g^-p on [g1, g2], gdot = -(a g + b g^2); run to t = 40/a (steady state)."""
    g = np.geomspace(2.0, 1e10, nb)
    Q = powerlaw(g, p, g1, g2)
    gdot = -(a * g + b * g ** 2)
    n = np.zeros(nb)
    for t_end, nsteps in ((1.0, 4000), (39.0, 4000)):   # fine steps first (fast cooling at high g), then coarse
        n = evolve(g, n, gdot, Q, t_end / nsteps, nsteps)
    cum = np.concatenate([np.cumsum((0.5 * (Q[1:] + Q[:-1]) * np.diff(g))[::-1])[::-1], [0.0]])  # int_g^inf Q
    exact = cum / np.abs(gdot)
    return g, n, exact


def rel_err(num, exact, sel):
    return np.max(np.abs(num[sel] / exact[sel] - 1.0))


NB_TEST = (550, 1100, 2200)   # numbins of the convergence study; 1100 = production default


def kernel_tests():
    """FP_FinDif_cool is first-order upwind in ln g: on a log grid a power law n ~ g^-s that is shifted in energy keeps
    its shape but its normalisation drifts at a relative rate ~ 1.3 dlng per e-fold of energy lost (s = 2.2;
    dlng = ln(g_max/g_min)/(numbins-1) = 0.020 for numbins = 1100), e.g. +8 % after the 3 e-folds of K2. The checks
    therefore bound the error at numbins = 1100 and require first-order convergence (error halves when numbins doubles)."""
    print("Kernel tests (FP_FinDif_cool)")
    out, conv = {}, {}
    cases = {
        "K1 impulsive, synchrotron": (lambda nb: k1_sync(nb), lambda g: (g > 2e3) & (g < 5e4), 0.01),
        "K2 impulsive, adiabatic (3 e-folds)": (lambda nb: k2_adiabatic(nb), lambda g: (g > 150) & (g < 1.6e6), 0.10),
        "K3 steady state, sync + adiabatic": (lambda nb: (lambda r: (r[0], None, r[1], r[2]))(k3_steady(nb)),
                                              lambda g: (g > 20) & (g < 3e7), 0.05),
    }
    for (name, (run, interior, tol)), key in zip(cases.items(), ("k1", "k2", "k3")):
        errs = {}
        for nb in NB_TEST:
            g, n0, n, exact = run(nb)
            errs[nb] = rel_err(n, exact, interior(g))
            if nb == 1100:
                out[key] = (g, n0, n, exact)
                if key == "k1":
                    check("K4 particle number conserved (synchrotron, no boundary loss)",
                          abs(number(g, n) / number(g, n0) - 1.0), 1e-10)
        conv[name] = errs
        order = np.log(errs[1100] / errs[2200]) / np.log(2200 / 1100)
        check(f"{name}: max rel. error in the interior, numbins=1100", errs[1100], tol)
        check(f"{name}: convergence order 1100->2200 >= 0.8 (shortfall)", max(0.0, 0.8 - order), 0.0)
        print(f"  info  {name}: " + ", ".join(f"{k}: {errs[k]:.2e}" for k in NB_TEST) + f"  (order {order:.2f})")
    out["k5"] = conv
    return out


# ----------------------------------------------------------------------------------------------------------------
# Full-model test
# ----------------------------------------------------------------------------------------------------------------
def scratch_namelist(src, outdir, numbins=None):
    txt = Path(src).read_text()
    if numbins is not None:
        txt, n = re.subn(r"\bnumbins\s*=\s*\d+", f"numbins = {numbins}", txt, flags=re.I)
        assert n == 1, "numbins not found in namelist"
    txt, n = re.subn(r"outdir\s*=\s*'[^']*'", f"outdir = '{outdir}'", txt)
    assert n == 1, "outdir not found in namelist"
    for key, val in (("CJ_", ".false."), ("radsyn_calc", ".false."), ("radpsyn_calc", ".false."),
                     ("radic_calc", ".false."), ("print_tlosses", ".false."), ("sort_outputs", ".false."),
                     ("NE_calc", ".true."), ("JET_", ".true.")):
        txt = re.sub(rf"\b{key}\s*=\s*\.\w+\.", f"{key} = {val}", txt, flags=re.I)
    path = Path(outdir) / "ne_validation.nml"
    path.write_text(txt)
    return path


class Model:
    """Loss and injection terms of the code, vectorised in g; spot-checked against the Fortran functions."""

    def __init__(self, nml):
        sp = J.system_parameters
        sp.read_namelist(str(nml))
        self.sp = sp
        self.z0, self.zmax, self.zinj = float(sp.z0), float(sp.zmax), float(sp.zmax_inj)
        self.Emin, self.Emax = float(sp.eminff), float(sp.emaxff)
        self.Emin_q = self.Emin if bool(J.controls.fix_inj_emin) else 1.1 * MEC2
        self.p = float(sp.p)
        self.gmin, self.gmax = self.Emin / MEC2, self.Emax / MEC2
        # per-z coefficients on a fine log-z table (interpolated in log-log)
        self.lz = np.linspace(np.log(self.z0), np.log(self.zmax * 1.01), 6000)
        zt = np.exp(self.lz)
        B = np.array([sp.bff_fun(z) for z in zt])
        bg = np.array([sp.betagamma(z) for z in zt])
        gam = np.array([sp.lorentzfactor(z) for z in zt])
        r2 = np.array([sp.r_fun_out(z) ** 2 - sp.r_fun_in(z) ** 2 for z in zt])
        beta = np.sqrt(1.0 - 1.0 / gam ** 2)
        # TLosses = SyncLoss + AdLoss  [erg/s]; SyncLoss = -1.05804e-15 B^2 (g^2-1), AdLoss = -(2/3) w Gamma beta c / z (E - mc2^2/E)
        self.cs = 1.05804e-15 * B ** 2 / MEC2 / (bg * C_LIGHT)                      # gdot_z sync  = -cs (g^2 - 1)
        self.ca = 2.0 / 3.0 * float(sp.geomindex_out) * gam * beta * C_LIGHT / zt / (bg * C_LIGHT)  # adiab = -ca (g - 1/g)
        Ecut = np.array([J.edistribution.cute(B[i], zt[i], 1.0, 1.0, "e") for i in range(len(zt))])
        q0p = np.array([J.edistribution.norm_q0p(Ecut[i], zt[i], "e") for i in range(len(zt))])
        self.lEcut = np.log(Ecut)
        self.lqlin = np.log(q0p * MEC2 * np.pi * r2) - zt / self.zinj      # Q_lin = Q_E * mc2 * pi r^2 * exp(-z/zinj)
        self.lbg, self.lr2 = np.log(bg), np.log(r2)
        self.lcs, self.lca = np.log(self.cs), np.log(self.ca)
        self._spot_check(zt, B, Ecut, q0p)

    def _i(self, arr, z):
        return np.interp(np.log(z), self.lz, arr)

    def gdot_z(self, g, z):
        return -(np.exp(self._i(self.lcs, z)) * (g ** 2 - 1.0) + np.exp(self._i(self.lca, z)) * (g - 1.0 / g))

    def dgdot_dg(self, g, z):
        return -(np.exp(self._i(self.lcs, z)) * 2.0 * g + np.exp(self._i(self.lca, z)) * (1.0 + 1.0 / g ** 2))

    def q_shape(self, E, Ecut):
        """Q_E / q0p (edistribution::Q_E with K = 1)."""
        x = E / Ecut
        q = (E / self.Emin_q) ** (-self.p) * np.where(x > 5e-4, np.exp(-x), 1.0)
        return np.where((E < self.Emin_q) | (E > self.Emax) | (x > 100.0), 0.0, q)

    def q_lin(self, g, z):
        return np.exp(self._i(self.lqlin, z)) * self.q_shape(g * MEC2, np.exp(self._i(self.lEcut, z)))

    def _spot_check(self, zt, B, Ecut, q0p):
        rng = np.random.default_rng(1)
        worst = 0.0
        for i in rng.integers(0, len(zt), 40):
            for E in np.geomspace(self.Emin * 1.01, self.Emax * 0.99, 25):
                g = E / MEC2
                ref = J.losses.tlosses(B[i], zt[i], 1.0, 1.0, E, "e") / MEC2 / (np.exp(self.lbg[i]) * C_LIGHT)
                mine = -(self.cs[i] * (g ** 2 - 1) + self.ca[i] * (g - 1 / g))
                worst = max(worst, abs(mine / ref - 1))
                qref = J.edistribution.q_e(q0p[i], Ecut[i], E, self.Emin_q, self.Emax)
                if qref > 0:
                    worst = max(worst, abs(q0p[i] * self.q_shape(E, Ecut[i]) / qref - 1))
        if worst > 1e-10:
            sys.exit(f"reference terms disagree with the Fortran TLosses/Q_E (max rel. diff {worst:.2e}): update Model")

    def to_ne(self, F, z):
        """F [per unit g] -> N_e per unit energy per unit comoving volume, as at the end of calcNE_paramo."""
        return F / (np.exp(self._i(self.lr2, z)) * np.pi * np.exp(self._i(self.lbg, z)) * C_LIGHT * MEC2)

    def characteristics(self, g_s, z_s):
        """F(g_s, z_s) for an array of g_s, integrating back along the cooling characteristics to z0.

        Per energy, in u = ln(z_s / z): d ln g'/du = -z gdot_z / g',  d ln J/du = -z dgdot_z/dg (J = dg'/dg),
        dF/du = z Q_lin(g', z) J. Backwards the energy only grows (not stiff); stop once g' leaves the injection range."""
        u_end = np.log(z_s / (self.z0 * 1.001))
        lgstop = np.log(self.gmax) + 0.01

        def rhs(u, y):
            z = z_s * np.exp(-u)
            g = np.exp(y[0])
            return [-z * self.gdot_z(g, z) / g, -z * self.dgdot_dg(g, z), z * self.q_lin(g, z) * np.exp(y[1])]

        def leave(u, y):
            return y[0] - lgstop
        leave.terminal, leave.direction = True, 1

        gq = np.geomspace(self.gmin, self.gmax, 2000)
        fscale = z_s * self.q_lin(gq, z_s).max()            # order of F for uncooled particles
        F = np.empty(len(g_s))
        for k, g0 in enumerate(g_s):
            sol = solve_ivp(rhs, (0.0, u_end), [np.log(g0), 0.0, 0.0], method="RK45", rtol=1e-7, atol=[1e-9, 1e-9, 1e-12 * fscale],
                            events=leave)
            assert sol.status >= 0, sol.message
            F[k] = sol.y[2, -1]
        return F

    def one_zone(self, g, z):
        gg = np.geomspace(self.gmin, self.gmax, 20000)
        q = self.q_lin(gg, z)
        cum = np.concatenate([np.cumsum((0.5 * (q[1:] + q[:-1]) * np.diff(gg))[::-1])[::-1], [0.0]])
        return np.interp(np.log(g), np.log(gg), cum) / np.abs(self.gdot_z(g, z))

    def g_cool(self, z):
        """g at which the synchrotron cooling length equals z (t_sync = t_adv = z / (betaGamma c))."""
        return 1.0 / (np.exp(self._i(self.lcs, z)) * z)


def run_ne(nml, outdir, numbins=None):
    """run.x with N_e only; returns (E [erg], N_e[E, segment], z of the segments, scratch namelist)."""
    outdir.mkdir(parents=True, exist_ok=True)
    snml = scratch_namelist(nml, outdir, numbins)
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    subprocess.run([str(CODE / "run.x"), str(snml)], cwd=CODE, env=env, check=True, stdout=subprocess.DEVNULL)
    ne = np.loadtxt(outdir / "Ne_cuts.dat")
    zseg = np.abs(np.loadtxt(outdir / "quant_FF.dat")[:, 0])
    E = ne[::len(zseg), 0]
    return E, ne[:, 1].reshape(len(E), len(zseg)), zseg, snml


def model_test(nml, heights_rg, workdir):
    """M1 at the namelist's numbins and at twice that. The kernel is first-order in ln g (see kernel_tests), so the
    error must halve; beyond zmax_inj (no injection, only cooling) it is largest, in the steep cooling cutoff, where
    numerical diffusion shifts the edge to higher energy. Bounded at 5 per cent only where injection is active."""
    nb0 = int(re.search(r"\bnumbins\s*=\s*(\d+)", Path(nml).read_text(), flags=re.I).group(1))
    print(f"Full-model test (run.x, N_e only) on {nml}, numbins = {nb0} and {2 * nb0}")
    E, N, zseg, snml = run_ne(nml, Path(workdir) / "run")
    _, N2, zseg2, _ = run_ne(nml, Path(workdir) / "run2", 2 * nb0)
    assert np.allclose(zseg, zseg2)
    M = Model(snml)
    rg = float(M.sp.r_g)
    g = E / MEC2
    res = []
    for h in heights_rg:
        j = int(np.argmin(np.abs(np.log(zseg[1:-1] / (h * rg))))) + 1      # segments 2..Ns-1 are filled
        z = zseg[j]
        ref = M.to_ne(M.characteristics(g, z), z)
        onez = M.to_ne(M.one_zone(g, z), z)
        w = E ** 2 * ref
        bulk = w > 1e-3 * w.max()                                         # where the energy is (E^2 N within 1e-3 of peak)
        fast = bulk & (g > 30 * M.g_cool(z))                              # t_sync < t_adv / 30
        e1, e2 = rel_err(N[:, j], ref, bulk), rel_err(N2[:, j], ref, bulk)
        order = np.log(e1 / e2) / np.log(2.0)
        res.append(dict(h=h, z=z, E=E, code=N[:, j], code2=N2[:, j], ref=ref, onez=onez, bulk=bulk, fast=fast,
                        gcool=M.g_cool(z)))
        tag = f"z = {z / rg:.3g} Rg" + ("" if z <= M.zinj else " (beyond zmax_inj)")
        print(f"  info  M1 {tag}: max |code/exact - 1| (E^2N > 1e-3 peak) = {e1:.2e} (numbins {nb0}), "
              f"{e2:.2e} ({2 * nb0}); order {order:.2f}")
        if z <= M.zinj:
            check(f"M1 N_e code vs characteristics, {tag}, numbins {nb0}", e1, 0.05)
        if e1 > 2e-3:
            check(f"M1 first-order convergence, {tag} (shortfall from 0.8)", max(0.0, 0.8 - order), 0.0)
        if fast.sum() > 3:
            print(f"  info  M2 one-zone vs characteristics where t_sync < t_adv/30: "
                  f"{rel_err(onez, ref, fast):.2e} ({fast.sum()} pts)")
    return M, res


# ----------------------------------------------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------------------------------------------
def plot_kernel(k, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 4, figsize=(16, 6.5), sharex="col", gridspec_kw=dict(height_ratios=(3, 1)))
    panels = [("K1  impulsive, synchrotron $\\dot\\gamma=-b\\gamma^2$", k["k1"][0], k["k1"][2], k["k1"][3], k["k1"][1]),
              ("K2  impulsive, adiabatic $\\dot\\gamma=-a\\gamma$", k["k2"][0], k["k2"][2], k["k2"][3], k["k2"][1]),
              ("K3  continuous injection, steady state", k["k3"][0], k["k3"][2], k["k3"][3], None)]
    for i, (title, g, n, ex, n0) in enumerate(panels):
        a = ax[0, i]
        if n0 is not None:
            a.loglog(g, g ** 2 * n0, color="0.7", lw=1, label="initial")
        a.loglog(g, g ** 2 * ex, "k-", lw=2.5, label="analytic")
        a.loglog(g, g ** 2 * n, "C3--", lw=1.5, label="FP_FinDif_cool")
        top = np.nanmax(g ** 2 * np.where(np.isfinite(ex), ex, 0))
        a.set_ylim(top * 1e-4, top * 3); a.set_title(title, fontsize=10); a.legend(fontsize=8)
        a.set_ylabel(r"$\gamma^2 n(\gamma)$")
        w = g ** 2 * np.where(np.isfinite(ex), ex, 0)
        s = w > 1e-3 * w.max()
        ax[1, i].semilogx(g[s], n[s] / ex[s] - 1, "C3")
        ax[1, i].axhline(0, color="k", lw=0.8); ax[1, i].set_ylim(-0.15, 0.15)
        ax[1, i].set_xlabel(r"$\gamma$"); ax[1, i].set_ylabel("num/exact − 1")
        a.set_xlim(2, 1e9)
    a = ax[0, 3]
    for j, (name, errs) in enumerate(k["k5"].items()):
        nbs = np.array(sorted(errs)); e = np.array([errs[n] for n in nbs])
        a.loglog(nbs, e, "o-", color=f"C{j}", label=name.split(":")[0])
        a.loglog(nbs, e[0] * nbs[0] / nbs, ":", color=f"C{j}", lw=0.8)
    a.plot([], [], "k:", lw=0.8, label=r"$\propto$ numbins$^{-1}$")
    for nb, lab in ((500, "paper, old"), (1100, "default")):
        a.axvline(nb, color="0.6", ls="--", lw=0.8)
        a.text(nb * 1.04, 0.3, lab, fontsize=8, rotation=90, transform=a.get_xaxis_transform())
    a.set_xticks([500, 1000, 2000]); a.set_xticklabels(["500", "1000", "2000"])
    a.set_xlabel("numbins"); a.set_ylabel("max |num/exact − 1| (interior)")
    a.set_title("K5  convergence with energy bins", fontsize=10); a.legend(fontsize=7)
    ax[1, 3].axis("off")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    print(f"  plot: {path}")


def plot_model(res, nml, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a, b) = plt.subplots(2, 1, figsize=(8.5, 8), sharex=True, gridspec_kw=dict(height_ratios=(3, 1.3)))
    eV = 1.602176634e-12
    np_err = np.seterr(divide="ignore", invalid="ignore")   # ratios where the reference is 0 are masked
    for i, r in enumerate(res):
        col = f"C{i}"
        E = r["E"] / eV                                   # x axis in eV; y = E^2 N with E in erg [erg cm^-3]
        Eerg = r["E"]
        lab = f"z ≈ {r['h']:.0e}".replace("e+0", "e") + r" $R_{\rm g}$"
        a.loglog(E, Eerg ** 2 * r["ref"], "-", color=col, lw=3, alpha=0.45)
        a.loglog(E, Eerg ** 2 * r["code"], "--", color=col, lw=1.4, label=lab)
        a.loglog(E, Eerg ** 2 * r["onez"], ":", color=col, lw=1.2)
        a.axvline(r["gcool"] * MEC2 / eV, color=col, lw=0.6, alpha=0.6)
        s = r["ref"] > 0
        b.semilogx(E[s & r["bulk"]], (r["code"] / r["ref"] - 1)[s & r["bulk"]], "-", color=col, lw=1.4)
        b.semilogx(E[s & r["bulk"]], (r["onez"] / r["ref"] - 1)[s & r["bulk"]], ":", color=col, lw=1.2)
        b.semilogx(E[s & r["bulk"]], (r["code2"] / r["ref"] - 1)[s & r["bulk"]], "-", color=col, lw=0.7, alpha=0.6)
    top = max((r["E"] ** 2 * r["ref"]).max() for r in res)
    a.set_ylim(top * 1e-7, top * 5)
    a.plot([], [], "k-", lw=3, alpha=0.45, label="characteristics (exact)")
    a.plot([], [], "k--", lw=1.4, label="code (calcNE_paramo)")
    a.plot([], [], "k:", lw=1.2, label="one-zone, local cooling")
    a.plot([], [], "k-", lw=0.7, alpha=0.6, label="code, 2× numbins (residuals)")
    a.plot([], [], "k-", lw=0.6, label=r"$t_{\rm sync}=t_{\rm adv}$")
    a.set_ylabel(r"$E^2 N_{\rm e}(E)$  [erg cm$^{-3}$]")
    a.set_title(f"N_e validation, {Path(nml).name}", fontsize=11)
    a.legend(fontsize=8, ncol=2)
    b.axhline(0, color="k", lw=0.8); b.set_ylim(-0.2, 0.6)
    b.set_xlabel(r"$E$ [eV] (comoving)"); b.set_ylabel("/ exact − 1")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    np.seterr(**np_err)
    print(f"  plot: {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nml", default=str(HERE.parent / "configs" / "m87_reference.nml"))
    ap.add_argument("--heights", type=float, nargs="+", default=[1e2, 1e3, 1e4, 1e5, 1e6], help="z in R_g")
    ap.add_argument("--plots", default=str(HERE / "output" / "ne_validation"))
    ap.add_argument("--skip-model", action="store_true")
    a = ap.parse_args()
    Path(a.plots).mkdir(parents=True, exist_ok=True)
    k = kernel_tests()
    plot_kernel(k, Path(a.plots) / "ne_kernel.png")
    if not a.skip_model:
        with tempfile.TemporaryDirectory(prefix="ne_validation_") as tmp:
            _, res = model_test(a.nml, a.heights, tmp)
        plot_model(res, a.nml, Path(a.plots) / f"ne_model_{Path(a.nml).stem}.png")
    failed = [r for r in RESULTS if not r[3]]
    print(f"\n{len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
