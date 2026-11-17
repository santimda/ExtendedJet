"""Bilby/emcee fitting runner."""

from __future__ import annotations

import shutil
from pathlib import Path

import bilby
import matplotlib.pyplot as plt
import numpy as np

from .likelihood import SEDLikelihood
from .parameters import build_priors
from .plotting import plot_sed


def fit_sed(
    model,
    data,
    fit_parameters,
    output_dir,
    fixed_parameters=None,
    truth=None,
    nwalkers=8,
    nsteps=200,
    seed=1,
    posterior_samples=200,
):
    fit_parameters = tuple(fit_parameters)
    if nwalkers < 2 * len(fit_parameters) or nsteps < 2:
        raise ValueError("use at least twice as many walkers as fitted parameters and two steps")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(output_dir / "chain", ignore_errors=True)
    priors = build_priors(fit_parameters)
    defaults = getattr(model, "parameters", {}).copy()
    defaults.update(fixed_parameters or {})
    likelihood = SEDLikelihood(model, data, fit_parameters, fixed_parameters=defaults)
    bilby.core.utils.random.seed(seed)
    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=priors,
        sampler="emcee",
        nwalkers=nwalkers,
        nsteps=nsteps,
        outdir=str(output_dir / "chain"),
        label="mcmc",
        npool=1,
        plot=False,
        resume=False,
    )
    corner = result.plot_corner(
        truth=None if truth is None else {
            name: (np.log10(truth["Lj"]) if name == "log_Lj" else truth[name])
            for name in fit_parameters
        },
        save=False,
    )
    corner.savefig(output_dir / "corner.pdf", bbox_inches="tight")
    plt.close(corner)
    plotting_seconds = plot_sed(
        likelihood,
        result.posterior,
        output_dir / "SED_MCMC.pdf",
        truth=truth,
        posterior_samples=posterior_samples,
    )
    print(f"MCMC model evaluations: {likelihood._evaluation} in {likelihood.evaluation_seconds:.1f}s")
    print(f"posterior plotting models: {posterior_samples + 1} in {plotting_seconds:.1f}s")
    print(f"wrote {output_dir / 'SED_MCMC.pdf'} and {output_dir / 'corner.pdf'}")
    return result
