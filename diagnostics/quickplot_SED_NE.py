#!/usr/bin/env python
"""Create the standard ExtendedJet particle and synchrotron diagnostic plots."""

import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from postprocessing.plot_utils import nice_fonts  # noqa: F401


OUTPUT_DIR = REPO_ROOT / "output"

eV = 1.60218e-12
h = 6.6260755e-27
mJy = 1e-26
pc = 3.086e18
D = 16.8e6 * pc
LUMINOSITY_FACTOR = 4 * np.pi * D**2


def save_figure(fig, filename):
    fig.savefig(filename, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def add_continuous_colorbar(fig, ax, cmap, norm, ticks, label):
    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(scalar_mappable, ticks=ticks, ax=ax, pad=0.01)
    colorbar.set_label(label)


def add_discrete_colorbar(fig, ax, cmap, count, label):
    bounds = np.arange(1, count + 2) - 0.5
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    scalar_mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(
        scalar_mappable,
        ticks=np.arange(1, count + 1),
        pad=0.01,
        ax=ax,
    )
    colorbar.set_label(label)


def plot_electron_cuts():
    """Plot the electron energy distribution for each jet segment."""
    data = np.loadtxt(OUTPUT_DIR / "Ne_cuts_ord.dat")
    energies = data[:, 0]
    densities = data[:, 1]
    segment_count = int(np.max(data[:, 2]))
    energy_count = len(data[:, 2]) // segment_count
    energies = energies.reshape(segment_count, energy_count)
    densities = densities.reshape(segment_count, energy_count)

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = matplotlib.colormaps["viridis"]
    norm = matplotlib.colors.Normalize(vmin=1, vmax=segment_count)
    for index in range(segment_count):
        ax.loglog(
            energies[index] / eV,
            densities[index] * energies[index] ** 2,
            color=cmap(norm(index + 1)),
        )

    add_continuous_colorbar(
        fig,
        ax,
        cmap,
        norm,
        np.linspace(0, segment_count, 11),
        "Segment index",
    )
    ax.set_xlim(1e6, 1e13)
    ax.set_ylim(1e-18, 1e-4)
    ax.set_xlabel(r"$E'$ [eV]")
    ax.set_ylabel(r"$E'^2 n'(E')$ [erg cm$^{-3}$]")
    save_figure(fig, "NEe_cuts.pdf")


def plot_total_electron_distribution():
    """Plot the electron distribution integrated up to each jet height."""
    data = np.loadtxt(OUTPUT_DIR / "Ne_bins_total.dat", skiprows=1)
    energies = data[:, 0]
    distributions = data[:, 1:-1]
    total = data[:, -1]
    bin_count = distributions.shape[1]

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = matplotlib.colormaps["viridis"].resampled(bin_count)
    colors = cmap(np.linspace(0, 1, bin_count))
    for index in range(bin_count):
        ax.loglog(
            energies,
            distributions[:, index] * (energies * eV) ** 2,
            color=colors[index],
            lw=1.5,
        )
    ax.loglog(
        energies,
        (energies * eV) ** 2 * total,
        color="darkred",
        lw=3,
        label="Total",
    )

    add_discrete_colorbar(fig, ax, cmap, bin_count, "Bin index")
    ax.set_xlim(1e6, 1e13)
    ax.set_ylim(1e42, 1e49)
    ax.set_xlabel(r"$E'$ [eV]")
    ax.set_ylabel(r"$E'^2 N'(E')$ [erg]")
    ax.legend(loc="best")
    save_figure(fig, "NEe_bins_total.pdf")


def read_synchrotron_cuts(filename):
    data = np.loadtxt(filename, skiprows=1).T
    energies, luminosities, segments = data
    segments = segments.astype(int)
    segment_ids = np.unique(segments)
    segment_count = len(segment_ids)
    energy_count = int(np.sum(segments == segment_ids[0]))
    energy_grid = energies[:energy_count]
    luminosity_grid = np.zeros((segment_count, energy_count))

    for index, segment in enumerate(segment_ids):
        mask = segments == segment
        luminosity_grid[index] = luminosities[mask] / energies[mask]

    return energy_grid, luminosity_grid, segment_count


def plot_synchrotron_cuts():
    """Plot the synchrotron SED for each jet segment."""
    energies, luminosities, segment_count = read_synchrotron_cuts(
        OUTPUT_DIR / "SEDe_syn_cuts_ord.dat"
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_xlim(2e-8, 1e9)
    ax.set_ylim(1e28, 1e39)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\epsilon$ [eV]")
    ax.set_ylabel(r"$\epsilon L_{\epsilon}$ [erg\,s$^{-1}$]")

    cmap = matplotlib.colormaps["viridis"]
    norm = matplotlib.colors.Normalize(vmin=1, vmax=segment_count)
    for index in range(segment_count):
        ax.loglog(
            energies,
            energies * luminosities[index],
            color=cmap(norm(index + 1)),
        )

    add_continuous_colorbar(
        fig,
        ax,
        cmap,
        norm,
        np.linspace(0, segment_count + 1, 11),
        "Segment index",
    )
    save_figure(fig, "SEDsyn_cuts.pdf")


def load_flux_outputs():
    jet = np.loadtxt(OUTPUT_DIR / "Flux_e.dat", skiprows=1)
    counter_jet = np.loadtxt(OUTPUT_DIR / "FluxCJ_e.dat", skiprows=1)
    return jet, counter_jet


def plot_luminosity_bins():
    """Plot the luminosity integrated up to each jet aperture."""
    jet, counter_jet = load_flux_outputs()
    energy = jet[:, 0] * h / eV
    counter_jet_energy = counter_jet[:, 0] * h / eV
    flux_bins = jet[:, 1:-1] / mJy
    flux_total = jet[:, -1] / mJy
    counter_jet_total = counter_jet[:, -1] / mJy
    bin_count = flux_bins.shape[1]

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = matplotlib.colormaps["viridis"].resampled(bin_count)
    colors = cmap(np.linspace(0, 1, bin_count))
    for index in range(bin_count):
        ax.loglog(
            energy,
            flux_bins[:, index] * mJy * jet[:, 0] * LUMINOSITY_FACTOR,
            color=colors[index],
            lw=1.5,
        )
    ax.loglog(
        energy,
        flux_total * mJy * jet[:, 0] * LUMINOSITY_FACTOR,
        color="darkred",
        lw=3,
        label="Total J",
    )
    ax.loglog(
        counter_jet_energy,
        counter_jet_total * mJy * counter_jet[:, 0] * LUMINOSITY_FACTOR,
        color="darkred",
        lw=3,
        linestyle="--",
        label="Total CJ",
    )

    add_discrete_colorbar(fig, ax, cmap, bin_count, "Bin index")
    ax.set_xlabel(r"$\epsilon$ [eV]")
    ax.set_ylabel(r"$\epsilon L_{\epsilon}$ [erg\,s$^{-1}$]")
    ax.set_xlim(1e-5, 1e7)
    ax.set_ylim(1e38, 1e41)
    ax.legend(loc="best")
    save_figure(fig, "Lum_e_bins.pdf")


def plot_radio_flux_bins():
    """Plot the flux density integrated up to each jet aperture."""
    jet, counter_jet = load_flux_outputs()
    frequency = jet[:, 0] / 1e9
    counter_jet_frequency = counter_jet[:, 0] / 1e9
    flux_bins = jet[:, 1:-1] / mJy
    flux_total = jet[:, -1] / mJy
    counter_jet_total = counter_jet[:, -1] / mJy
    bin_count = flux_bins.shape[1]

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = matplotlib.colormaps["viridis"].resampled(bin_count)
    colors = cmap(np.linspace(0, 1, bin_count))
    for index in range(bin_count):
        ax.loglog(frequency, flux_bins[:, index], color=colors[index], lw=1.5)
    ax.loglog(frequency, flux_total, color="darkred", lw=3, label="Total J")
    ax.loglog(
        counter_jet_frequency,
        counter_jet_total,
        color="darkred",
        lw=3,
        linestyle="--",
        label="Total CJ",
    )

    add_discrete_colorbar(fig, ax, cmap, bin_count, "Bin index")
    ax.set_xlabel(r"$\nu$ [GHz]")
    ax.set_ylabel(r"$S_{\nu}$ [mJy]")
    ax.set_xlim(1e-2, 1e3)
    ax.set_ylim(1e-1, 5e3)
    ax.legend(loc="best")
    save_figure(fig, "Flux_e_bins_radio.pdf")


def main():
    plot_electron_cuts()
    plot_total_electron_distribution()
    plot_synchrotron_cuts()
    plot_luminosity_bins()
    plot_radio_flux_bins()


if __name__ == "__main__":
    main()
