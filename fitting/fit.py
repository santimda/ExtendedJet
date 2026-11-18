"""Common priors, likelihood, sampling, and output generation."""

from __future__ import annotations

import shutil
import time
from pathlib import Path

import bilby
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from postprocessing.plot_utils import nice_fonts

from .config import FluxData
from .plotting import plot_sed

matplotlib.rcParams.update(nice_fonts)

PRIOR_DEFINITIONS = {
    "log_Lj": {
        "bounds": (40.0, 45.0),
        "latex_label": r"$\log_{10}(L_j/[\mathrm{erg\ s^{-1}}])$",
    },
    "log_epse": {
        "bounds": (-3.0, -0.5),
        "latex_label": r"$\log_{10}(\epsilon_e)$",
    },
    "log_z0": {
        "bounds": (float(np.log10(2.0)), float(np.log10(200.0))),
        "latex_label": r"$\log_{10}(z_0/R_g)$",
    },
    "log_zinj": {
        "bounds": (2.5, 5.5),
        "latex_label": r"$\log_{10}(z_{\rm inj}/R_g)$",
    },
    "g": {"bounds": (0.1, 0.5), "latex_label": r"$g$"},
    "log_sig0": {
        "bounds": (float(np.log10(2.0)), 3.0),
        "latex_label": r"$\log_{10}(\sigma_0)$",
    },
    "log_etaacc": {
        "bounds": (-4.0, -1.0),
        "latex_label": r"$\log_{10}(\eta_{\rm acc})$",
    },
    "p": {"bounds": (1.9, 3.0), "latex_label": r"$p$"},
    "geomIndex_out": {"bounds": (0.0, 1.0), "latex_label": r"$w$"},
    "a_inj": {"bounds": (-0.5, 1.0), "latex_label": r"$a_e$"},
}


def build_priors(fit_parameters):
    fit_parameters = tuple(fit_parameters)
    unknown = set(fit_parameters) - set(PRIOR_DEFINITIONS)
    if unknown:
        raise ValueError(f"No common prior is defined for: {sorted(unknown)}")
    priors = bilby.core.prior.PriorDict()
    for name in fit_parameters:
        definition = PRIOR_DEFINITIONS[name]
        minimum, maximum = definition["bounds"]
        priors[name] = bilby.core.prior.Uniform(
            minimum,
            maximum,
            name=name,
            latex_label=definition["latex_label"],
        )
    return priors


class SEDLikelihood(bilby.Likelihood):
    def __init__(self, model, data: FluxData, fixed_parameters=None):
        super().__init__()
        self.model = model
        self.data = data
        self.frequencies = data.frequencies_ghz
        self.apertures = data.apertures_rg
        self.fit_parameters = data.fit_parameters
        self.fixed_parameters = model.parameters.copy()
        self.fixed_parameters.update(fixed_parameters or {})
        self.observed = data.observed
        self.sigma = data.errors_cgs
        self._evaluation = 0
        self.evaluation_seconds = 0.0

    def evaluate_model(self, parameters, frequencies=None, apertures=None, output_dir=None):
        model_parameters = self.fixed_parameters.copy()
        model_parameters.update(parameters)
        for public_name, internal_name in self.model.log_parameters.items():
            if public_name in model_parameters:
                model_parameters.pop(internal_name, None)
        for public_name, internal_name in self.model.parameter_aliases.items():
            if public_name in model_parameters:
                model_parameters.pop(internal_name, None)
        return self.model(
            self.frequencies if frequencies is None else frequencies,
            self.apertures if apertures is None else apertures,
            model_parameters,
            output_dir=output_dir,
        )

    def log_likelihood(self, parameters):
        self._evaluation += 1
        started = time.perf_counter()
        try:
            prediction = self.evaluate_model(parameters)
        except (OSError, RuntimeError, ValueError):
            return -1e30
        finally:
            self.evaluation_seconds += time.perf_counter() - started
        residual = (self.data.fluxes_cgs - prediction)[self.observed] / self.sigma[self.observed]
        sigma = self.sigma[self.observed]
        return float(-0.5 * np.sum(residual**2 + np.log(2.0 * np.pi * sigma**2)))


def fit_sed(model, data, output_dir, mode, fixed_parameters=None, truth=None,
            nwalkers=8, nsteps=200, seed=1, posterior_samples=100):
    if mode not in {"overwrite", "resume"}:
        raise ValueError("mode must be 'overwrite' or 'resume'")
    if nwalkers < 2 * len(data.fit_parameters) or nsteps < 2:
        raise ValueError("use at least twice as many walkers as fitted parameters and two steps")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    chain_dir = output_dir / "chain"
    if mode == "overwrite":
        shutil.rmtree(chain_dir, ignore_errors=True)
    resume = mode == "resume" and (chain_dir / "emcee_mcmc" / "sampler.pickle").exists()
    if mode == "resume" and not resume:
        print("No existing emcee chain found; initializing a new chain.")
    priors = build_priors(data.fit_parameters)
    likelihood = SEDLikelihood(model, data, fixed_parameters=fixed_parameters)
    bilby.core.utils.random.seed(seed)
    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=priors,
        sampler="emcee",
        nwalkers=nwalkers,
        nsteps=nsteps,
        outdir=str(chain_dir),
        label="mcmc",
        npool=1,
        plot=False,
        resume=resume,
    )
    matplotlib.rcParams.update(nice_fonts)
    corner = result.plot_corner(
        truth=None if truth is None else {
            name: truth[name]
            for name in data.fit_parameters
        },
        save=False,
    )
    source = model.source_name
    corner.savefig(output_dir / f"corner_{source}.pdf", bbox_inches="tight")
    plt.close(corner)
    plotting_seconds = plot_sed(
        likelihood,
        result.posterior,
        output_dir / f"SED_{source}.pdf",
        truth=truth,
        posterior_samples=posterior_samples,
    )
    sample_count = min(posterior_samples, len(result.posterior))
    print(f"MCMC model evaluations: {likelihood._evaluation} in {likelihood.evaluation_seconds:.1f}s")
    print(f"posterior plotting models: {sample_count + 2 + (truth is not None)} in {plotting_seconds:.1f}s")
    print(f"wrote {output_dir / f'SED_{source}.pdf'} and {output_dir / f'corner_{source}.pdf'}")
    return result
