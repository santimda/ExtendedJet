"""Minimal Bilby/emcee recovery test for two ExtendedJet parameters.

The default run generates noiseless synthetic fluxes at five frequencies from
3 to 86 GHz and at apertures 30, 100, and 200 R_g. It then fits ``a_inj``
(``a_e``) and ``geomIndex_out`` (``w``) with a reduced chain.

    python validation/mcmc_sanity.py
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import bilby
import matplotlib.pyplot as plt
import numpy as np

try:
    from postprocessing import plot_utils  # noqa: F401
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from postprocessing import plot_utils  # noqa: F401


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NML = ROOT / "configs" / "m87_reference.nml"
DEFAULT_OUTDIR = ROOT / "validation" / "output" / "mcmc_validation"
FREQUENCIES_GHZ = np.array([3.0, 10.0, 22.0, 43.0, 86.0])
APERTURES_RG = np.array([30.0, 100.0, 200.0])
ERROR_FRACTION = 0.15
PARAMETERS = ("a_inj", "geomIndex_out")


def _replace(text, key, value):
    pattern = rf"(^\s*{re.escape(key)}\s*=\s*)[^,\n]+"
    updated, count = re.subn(pattern, rf"\g<1>{value}", text, flags=re.IGNORECASE | re.MULTILINE)
    if count != 1:
        raise ValueError(f"Expected one namelist value for {key}, found {count}")
    return updated


def make_namelist(base, outdir, parameters, ns, numbins, apertures):
    text = Path(base).read_text()
    for key, value in parameters.items():
        text = _replace(text, key, f"{value:.12g}d0")
    text = _replace(text, "outdir", f"'{Path(outdir).resolve()}'")
    text = _replace(text, "Ns", str(ns))
    text = _replace(text, "numbins", str(numbins))
    aperture_values = ", ".join(f"{value:.12g}d0" for value in apertures)
    aperture_values += ", " + ", ".join("0d0" for _ in range(10 - len(apertures)))
    text, count = re.subn(r"^\s*apertures_rg\s*=.*$", f"  apertures_rg = {aperture_values}", text, flags=re.IGNORECASE | re.MULTILINE)
    if count != 1:
        raise ValueError("Expected one apertures_rg entry")
    text = _replace(text, "n_ap", str(len(apertures)))
    text = _replace(text, "sort_outputs", ".false.")
    text = _replace(text, "print_tlosses", ".false.")
    path = Path(outdir) / "run.nml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def run_model(base, workdir, parameters, ns, numbins, apertures):
    workdir = Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    nml = make_namelist(base, workdir, parameters, ns, numbins, apertures)
    run_x = ROOT / "Code" / "run.x"
    process = subprocess.run([str(run_x), str(nml)], cwd=ROOT, capture_output=True, text=True)
    if process.returncode != 0:
        raise RuntimeError(f"run.x failed for {parameters}:\n{process.stdout}\n{process.stderr}")
    return load_flux(workdir)


def load_flux(workdir):
    jet = np.loadtxt(Path(workdir) / "Flux_e.dat", comments="#")
    counterjet = np.loadtxt(Path(workdir) / "FluxCJ_e.dat", comments="#")
    if jet.shape != counterjet.shape or jet.shape[1] < 4:
        raise ValueError("Jet and counter-jet flux tables are incompatible")
    return jet[:, 0], jet[:, 1:4] + counterjet[:, 1:4]


def interpolate_flux(nu_hz, flux, frequencies_ghz):
    requested = np.asarray(frequencies_ghz, float) * 1e9
    if requested[0] < nu_hz[0] or requested[-1] > nu_hz[-1]:
        raise ValueError("Requested frequencies are outside the model SED grid")
    curves = []
    for index in range(flux.shape[1]):
        positive = flux[:, index] > 0
        if not np.any(positive):
            curves.append(np.zeros_like(requested))
            continue
        positive_nu = nu_hz[positive]
        positive_flux = flux[positive, index]
        curves.append(np.where(
            requested < positive_nu[0],
            0.0,
            np.exp(np.interp(np.log(requested), np.log(positive_nu), np.log(positive_flux))),
        ))
    return np.array(curves).T


class SyntheticSEDLikelihood(bilby.Likelihood):
    def __init__(self, base, workdir, truth, frequencies, apertures, ns, numbins):
        super().__init__()
        self.base = Path(base)
        self.workdir = Path(workdir)
        self.truth = truth
        self.frequencies = np.asarray(frequencies, float)
        self.apertures = np.asarray(apertures, float)
        self.ns = ns
        self.numbins = numbins
        self.nu_hz, self.truth_curve = run_model(self.base, self.workdir / "truth", truth, ns, numbins, self.apertures)
        self.data = interpolate_flux(self.nu_hz, self.truth_curve, self.frequencies)
        self.sigma = np.maximum(ERROR_FRACTION * self.data, 1e-6 * np.max(self.data))
        self._evaluation = 0
        self.evaluation_seconds = 0.0

    def log_likelihood(self, parameters):
        model_parameters = {name: float(parameters[name]) for name in PARAMETERS}
        evaluation = self._evaluation
        self._evaluation += 1
        started = time.perf_counter()
        try:
            nu_hz, curve = run_model(self.base, self.workdir / f"eval_{evaluation:04d}", model_parameters, self.ns, self.numbins, self.apertures)
            prediction = interpolate_flux(nu_hz, curve, self.frequencies)
        except (OSError, RuntimeError, ValueError):
            return -1e30
        finally:
            self.evaluation_seconds += time.perf_counter() - started
        residual = (self.data - prediction) / self.sigma
        return float(-0.5 * np.sum(residual ** 2 + np.log(2.0 * np.pi * self.sigma ** 2)))


def plot_sed(likelihood, posterior, output):
    best = {name: float(posterior[name].median()) for name in PARAMETERS}
    dense_frequencies = np.geomspace(1.0, 100.0, 200)
    started = time.perf_counter()
    nu_hz, fitted_curve = run_model(likelihood.base, likelihood.workdir / "best", best, likelihood.ns, likelihood.numbins, likelihood.apertures)
    fitted = interpolate_flux(nu_hz, fitted_curve, dense_frequencies)
    truth = interpolate_flux(likelihood.nu_hz, likelihood.truth_curve, dense_frequencies)

    posterior_samples = posterior[list(PARAMETERS)].sample(n=200, random_state=1)
    sample_curves = []
    for index, sample in enumerate(posterior_samples.to_dict("records")):
        nu_hz, curve = run_model(
            likelihood.base,
            likelihood.workdir / f"posterior_{index:03d}",
            sample,
            likelihood.ns,
            likelihood.numbins,
            likelihood.apertures,
        )
        sample_curves.append(interpolate_flux(nu_hz, curve, dense_frequencies))
    sample_curves = np.asarray(sample_curves)
    lower = np.percentile(sample_curves, 16, axis=0)
    upper = np.percentile(sample_curves, 84, axis=0)

    frequencies = likelihood.frequencies
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(likelihood.apertures)))
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": (3, 1), "hspace": 0.03})
    for index, (aperture, color) in enumerate(zip(likelihood.apertures, colors)):
        data = likelihood.data[:, index]
        error = likelihood.sigma[:, index]
        model = interpolate_flux(nu_hz, fitted_curve, frequencies)[:, index]
        axes[0].plot(dense_frequencies, truth[:, index] / 1e-23, color=color, ls="--")
        axes[0].plot(dense_frequencies, fitted[:, index] / 1e-23, color=color)
        axes[0].fill_between(dense_frequencies, lower[:, index] / 1e-23, upper[:, index] / 1e-23, color=color, alpha=0.25, edgecolor="none")
        axes[0].errorbar(
            frequencies,
            data / 1e-23,
            yerr=error / 1e-23,
            fmt="o",
            markerfacecolor="none",
            color=color,
            alpha=1.0,
            markersize=6,
            elinewidth=1,
            zorder=12,
            label=rf"{aperture:g} $R_g$",
        )
        density_residual = np.divide(
            data - model, model, out=np.zeros_like(model), where=model != 0,
        )
        axes[1].plot(frequencies, density_residual, color=color, marker="o", ls="")
    axes[0].set(xscale="log", yscale="log", ylabel=r"$S_\nu$ [Jy]")
    axes[1].set(xscale="log", xlabel="Frequency [GHz]", ylabel=r"$\Delta S/S$", ylim=(-1.0, 1.0))
    axes[0].legend(fontsize=10)
    axes[1].axhline(0.0, color="k", ls="--", lw=1)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return time.perf_counter() - started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nml", type=Path, default=DEFAULT_NML)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--frequencies", type=float, nargs="+", default=FREQUENCIES_GHZ)
    parser.add_argument("--apertures", type=float, nargs="+", default=APERTURES_RG)
    parser.add_argument("--ns", type=int, default=40)
    parser.add_argument("--numbins", type=int, default=100)
    parser.add_argument("--nwalkers", type=int, default=8)
    parser.add_argument("--nsteps", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    if len(args.apertures) != 3 or np.any(np.asarray(args.apertures) <= 0):
        parser.error("exactly three positive apertures are required")
    if args.nwalkers < 4 or args.nsteps < 2:
        parser.error("use at least four walkers and two steps")
    args.outdir.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(args.outdir / "chain", ignore_errors=True)
    truth = {"a_inj": 0.01, "geomIndex_out": 0.56}
    likelihood = SyntheticSEDLikelihood(args.nml, args.outdir / "runs", truth, args.frequencies, args.apertures, args.ns, args.numbins)
    priors = bilby.core.prior.PriorDict({
        "a_inj": bilby.core.prior.Uniform(-0.5, 2, name="a_inj", latex_label=r"$a_e$"),
        "geomIndex_out": bilby.core.prior.Uniform(0.0, 1.0, name="geomIndex_out", latex_label=r"$w$"),
    })
    bilby.core.utils.random.seed(args.seed)
    result = bilby.run_sampler(likelihood=likelihood, priors=priors, sampler="emcee", nwalkers=args.nwalkers, nsteps=args.nsteps, outdir=str(args.outdir / "chain"), label="mcmc", npool=1, plot=False, resume=False)
    posterior = result.posterior
    corner_figure = result.plot_corner(
        truth={"a_inj": truth["a_inj"], "geomIndex_out": truth["geomIndex_out"]},
        save=False,
    )
    corner_figure.savefig(args.outdir / "corner.png", dpi=160, bbox_inches="tight")
    plt.close(corner_figure)
    plotting_seconds = plot_sed(likelihood, posterior, args.outdir / "SED_MCMC.pdf")
    print(f"MCMC model evaluations: {likelihood._evaluation} in {likelihood.evaluation_seconds:.1f}s")
    print(f"posterior plotting models: 201 in {plotting_seconds:.1f}s")
    print(f"wrote {args.outdir / 'SED_MCMC.pdf'} and {args.outdir / 'corner.png'}")


if __name__ == "__main__":
    main()
