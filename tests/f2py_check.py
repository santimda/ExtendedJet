#!/usr/bin/env python3
"""
Regression check: does the f2py `syn_model` wrapper (JetCode module, built via
`make f2py`) reproduce the compiled `run.x` binary to better than 1% relative
difference in flux, over 1e9-1e20 Hz?

Procedure
---------
1. Run the compiled binary `run.x` on the paper's Jet+RIAF namelist
   (analysis/configs/jet_riaf.nml, only `outdir` redirected to a scratch dir)
   and read the reference total flux (jet + counter-jet, no aperture cut) from
   `Flux_e.dat` and `FluxCJ_e.dat` (last column of each, summed).
2. Call `JetCode.system_parameters.read_namelist()` on the SAME namelist file,
   so every module-level parameter (including the `fix_norm_taper` /
   `fix_inj_emin` controls in module `controls`, which default to `.false.`
   and are NOT arguments of `syn_model`) matches exactly what `run.x` used.
3. Call `JetCode.main.syn_model(...)` with the derived (cgs) values that
   `read_namelist` just computed, at photon energies E_data = h*nu on the
   run.x frequency grid restricted to [1e9, 1e20] Hz.
4. Unit conversion: `synrad.f90` shows `SEDsyn` (which fills Flux_e.dat)
   accumulates `Total_flux = Total_flux + dFdz*h`, while `SEDsyn_data` (used
   by syn_model) accumulates `S_data = S_data + dFdz` -- i.e. WITHOUT the `h`
   factor. So Flux_e.dat's tabulated flux == h * syn_model's S_data at the
   same E_data. We multiply by h before comparing.
5. Report the max relative difference and the frequency at which it occurs,
   and the wall-clock time of a single syn_model call. Exit 1 if max relative
   difference > 1%.

Run from anywhere with the m87 conda env active:
    python codes/extended-jet/tests/f2py_check.py
"""
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(HERE, "..", "Code"))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
NML_SRC = os.path.join(REPO_ROOT, "analysis", "configs", "jet_riaf.nml")

FREQ_LO, FREQ_HI = 1e9, 1e20
REL_TOL = 0.01


def make_scratch_nml(scratch_dir):
    """Copy the paper's Jet+RIAF namelist, redirecting outdir to scratch_dir/outdir."""
    outdir = os.path.join(scratch_dir, "outdir")
    os.makedirs(outdir, exist_ok=True)
    with open(NML_SRC) as f:
        text = f.read()
    nml_path = os.path.join(scratch_dir, "jet_riaf_ref.nml")
    replaced = False
    lines = []
    for line in text.splitlines():
        if "outdir" in line and "=" in line and not replaced:
            lines.append(f"  outdir = '{outdir}'")
            replaced = True
        else:
            lines.append(line)
    if not replaced:
        raise RuntimeError("Could not find outdir line in namelist to replace")
    with open(nml_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return nml_path, outdir


def run_reference(nml_path):
    """Run the compiled run.x binary and return (nu[Hz], total_flux) on its full grid."""
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    subprocess.run(
        ["./run.x", nml_path], cwd=CODE_DIR, env=env, check=True,
        capture_output=True, text=True,
    )
    outdir = None
    with open(nml_path) as f:
        for line in f:
            if "outdir" in line and "=" in line:
                outdir = line.split("=", 1)[1].strip().strip(",").strip("'\"")
                break
    flux_j = np.loadtxt(os.path.join(outdir, "Flux_e.dat"))
    flux_cj = np.loadtxt(os.path.join(outdir, "FluxCJ_e.dat"))
    nu = flux_j[:, 0]
    assert np.allclose(nu, flux_cj[:, 0]), "J and CJ frequency grids differ"
    total = flux_j[:, -1] + flux_cj[:, -1]  # last column = no aperture cut, full jet
    return nu, total


def run_f2py(nml_path, nu):
    """Call JetCode.syn_model with parameters read from the same namelist."""
    sys.path.insert(0, CODE_DIR)
    import JetCode  # noqa: E402  (built by `make f2py`)

    sp = JetCode.system_parameters
    ct = JetCode.controls

    sp.read_namelist(nml_path)

    # Document the module state that syn_model will use, so it can be checked
    # against the &model / &run groups in the namelist used for the reference run.
    print("=== module state after read_namelist() (must match nml used for run.x) ===")
    print(f"sigma_0       = {sp.sigma_0!r}")
    print(f"zmax [cm]     = {sp.zmax!r}  (zmax_rg={sp.zmax_rg!r}, R_g={sp.r_g!r})")
    print(f"Lj            = {sp.lj!r}")
    print(f"a_inj         = {sp.a_inj!r}")
    print(f"zmax_inj [cm] = {sp.zmax_inj!r}  (zmax_inj_rg={sp.zmax_inj_rg!r})")
    print(f"gammaIndex    = {sp.gammaindex!r}")
    print(f"geomIndex_out = {sp.geomindex_out!r}")
    print(f"z0 [cm]       = {sp.z0!r}  (z0_rg={sp.z0_rg!r})")
    print(f"fnT           = {sp.fnt!r}")
    print(f"facc          = {sp.facc!r}")
    print(f"p             = {sp.p!r}")
    print(f"fix_norm_taper (controls) = {ct.fix_norm_taper!r}")
    print(f"fix_inj_emin   (controls) = {ct.fix_inj_emin!r}")
    print("=" * 78)

    h = JetCode.constants.h  # Planck's constant, cgs -- same value used in synrad.f90
    E_data = h * nu  # nu -> photon energy [erg]

    t0 = time.perf_counter()
    s_data = JetCode.syn_model(
        sp.sigma_0, sp.zmax, sp.lj, sp.a_inj,
        sp.zmax_inj, sp.gammaindex, sp.geomindex_out,
        sp.z0, sp.fnt, sp.facc, sp.p,
        E_data,
    )
    elapsed = time.perf_counter() - t0

    n = len(nu)
    s_data = np.asarray(s_data)[:n]  # only the first len(E_data) entries are filled

    # Flux_e.dat's Total_flux = h * S_data (see module docstring above).
    flux_f2py = h * s_data
    return flux_f2py, elapsed


def main():
    scratch_dir = os.environ.get("F2PY_CHECK_SCRATCH") or tempfile.mkdtemp(prefix="f2py_check_")
    os.makedirs(scratch_dir, exist_ok=True)
    nml_path, _ = make_scratch_nml(scratch_dir)

    print(f"Reference run: ./run.x {nml_path}")
    nu_full, flux_ref_full = run_reference(nml_path)

    mask = (nu_full >= FREQ_LO) & (nu_full <= FREQ_HI)
    nu = nu_full[mask]
    flux_ref = flux_ref_full[mask]
    print(f"Comparing on {len(nu)} points in [{FREQ_LO:.1e}, {FREQ_HI:.1e}] Hz "
          f"(of {len(nu_full)} total grid points)")

    flux_f2py_full, elapsed = run_f2py(nml_path, nu_full)
    flux_f2py = flux_f2py_full[mask]

    print(f"syn_model wall-clock time (single call, {len(nu_full)} energies): "
          f"{elapsed:.3f} s")

    # Relative difference, guarding against zero reference values.
    nonzero = flux_ref != 0
    reldiff = np.full_like(flux_ref, np.nan)
    reldiff[nonzero] = np.abs(flux_f2py[nonzero] - flux_ref[nonzero]) / np.abs(flux_ref[nonzero])
    if np.any(~nonzero):
        # where reference is exactly zero, flag any nonzero f2py output as a hard mismatch
        reldiff[~nonzero] = np.where(flux_f2py[~nonzero] == 0, 0.0, np.inf)

    imax = np.nanargmax(reldiff)
    max_reldiff = reldiff[imax]
    freq_at_max = nu[imax]

    print(f"Max relative difference: {max_reldiff:.6e} at nu = {freq_at_max:.6e} Hz")
    print(f"  run.x  flux at that point: {flux_ref[imax]:.6e}")
    print(f"  f2py   flux at that point: {flux_f2py[imax]:.6e}")

    if max_reldiff > REL_TOL:
        print(f"FAIL: max relative difference {max_reldiff:.4%} exceeds tolerance {REL_TOL:.0%}")
        sys.exit(1)
    else:
        print(f"PASS: max relative difference {max_reldiff:.4%} within tolerance {REL_TOL:.0%}")
        sys.exit(0)


if __name__ == "__main__":
    main()
