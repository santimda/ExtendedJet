"""SED and posterior plotting."""

from __future__ import annotations

import time

import matplotlib.pyplot as plt
import numpy as np

from .config import Jy


def plot_sed(likelihood, posterior, output, truth=None, posterior_samples=100):
    best = {name: float(posterior[name].median()) for name in likelihood.fit_parameters}
    dense_frequencies = np.geomspace(1.0, 1000.0, 200)
    started = time.perf_counter()
    fitted = likelihood.evaluate_model(best, frequencies=dense_frequencies, output_dir=likelihood.model.output_dir)
    fitted_data = likelihood.evaluate_model(best, output_dir=likelihood.model.output_dir)
    truth_curve = None
    if truth is not None:
        truth_curve = likelihood.evaluate_model(truth, frequencies=dense_frequencies, output_dir=likelihood.model.output_dir)

    samples = posterior[list(likelihood.fit_parameters)].sample(
        n=min(posterior_samples, len(posterior)), random_state=1
    )
    sample_curves = np.asarray([
        likelihood.evaluate_model(
            sample,
            frequencies=dense_frequencies,
            output_dir=likelihood.model.output_dir,
        )
        for sample in samples.to_dict("records")
    ])
    lower = np.percentile(sample_curves, 16, axis=0)
    upper = np.percentile(sample_curves, 84, axis=0)

    frequencies = likelihood.frequencies
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(likelihood.apertures)))
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8, 7),
        sharex=True,
        gridspec_kw={"height_ratios": (3, 0.7), "hspace": 0.03},
    )
    for index, (aperture, color) in enumerate(zip(likelihood.apertures, colors)):
        observed = likelihood.observed[:, index]
        data = likelihood.data.fluxes_cgs[:, index]
        error = likelihood.sigma[:, index]
        model = fitted_data[:, index]
        axes[0].plot(dense_frequencies, fitted[:, index] / Jy, color=color)
        if truth_curve is not None:
            axes[0].plot(dense_frequencies, truth_curve[:, index] / Jy, color=color, ls="--")
        axes[0].fill_between(
            dense_frequencies,
            lower[:, index] / Jy,
            upper[:, index] / Jy,
            color=color,
            alpha=0.25,
            edgecolor="none",
        )
        axes[0].errorbar(
            frequencies[observed],
            data[observed] / Jy,
            yerr=error[observed] / Jy,
            fmt="o",
            markerfacecolor="none",
            color=color,
            markersize=6,
            elinewidth=1,
            zorder=12,
            label=rf"{aperture:g} $R_g$",
        )
        residual = np.divide(data - model, model, out=np.full_like(model, np.nan), where=model != 0)
        axes[1].plot(frequencies[observed], residual[observed], color=color, marker="o", ls="")
    axes[0].set(xscale="log", yscale="log", xlim=(1.0, 1000.0), ylim=(0.01, 10.0), ylabel=r"$S_\nu$ [Jy]")
    axes[1].set(xscale="log", xlabel="Frequency [GHz]", ylabel=r"$\Delta S/S$", ylim=(-1.0, 1.0))
    axes[0].legend(fontsize=10)
    axes[1].axhline(0.0, color="k", ls="--", lw=1)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return time.perf_counter() - started
