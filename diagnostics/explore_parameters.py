#!/usr/bin/env python3
"""Explore ExtendedJet parameters with the f2py model wrapper.

Build the wrapper first with ``make -C Code f2py``, then run, for example::

    OMP_NUM_THREADS=1 python diagnostics/explore_parameters.py

The default reference values and sweep grids intentionally match the old
``explore_parameters.ipynb`` notebook.  Model calls are sequential because
the Fortran modules contain mutable global state.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from postprocessing.plot_utils import nice_fonts

matplotlib.rcParams.update(nice_fonts)


REPO_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = REPO_ROOT / "Code"
sys.path.insert(0, str(CODE_DIR))

eV = 1.602176e-12
h = 6.6260755e-27
M_BH = 6.5e9
R_g = 1.477e5 * M_BH
DIST = 16.8 * 3.0856e24
LUMINOSITY_FACTOR = 4.0 * np.pi * DIST**2

REFERENCE = {
    "sigma0": 1.0e2,
    "zmax": 1.1e5 * R_g,
    "Lj": 1.0e43,
    "a_inj": 0.01,
    "zmax_inj": 1.0e4 * R_g,
    "gammaIndex": 0.1,
    "geomIndex_out": 0.56,
    "z0": 10.0 * R_g,
    "fnT": 0.1,
    "facc": 1.0e-2,
    "p": 2.05,
}

E = np.geomspace(1.0e-7 * eV, 6.0e7 * eV, 80)
Z = np.geomspace(REFERENCE["z0"], REFERENCE["zmax"], 80)
SED_XMIN, SED_XMAX = E[0] / h, E[-1] / h
SED_YMIN, SED_YMAX = 1.0e39, 1.0e42
RADIO_XMIN, RADIO_XMAX = 1.0e-1, 1.0e3
RADIO_YMIN, RADIO_YMAX = 5.0e-2, 1.0e2

SWEEPS = {
    "a_inj": {
        "values": [-0.5, 0.01, 0.5, 0.99],
        "label": r"$a_{\rm inj}$",
        "filename": "SEDs_a_injs.pdf",
        "format": ".2g",
    },
    "zmax": {
        "values": np.array([1.2e5, 2e5, 3e5, 5e5, 1e7]) * R_g,
        "label": r"$z_{\rm max}/R_g$",
        "filename": "SEDs_zmax.pdf",
        "format": ".1e",
    },
    "sigma0": {
        "values": [2.0, 5.0, 30.0, 100.0],
        "label": r"$\sigma_0$",
        "filename": "SEDs_sigma0s.pdf",
        "format": ".2g",
    },
    "Lj": {
        "values": [1e42, 3e42, 1e43, 3e43],
        "label": r"$L_{\rm j}$ [erg s$^{-1}$]",
        "filename": "SEDs_Ljs.pdf",
        "format": ".1e",
    },
    "z0": {
        "values": np.array([2.0, 10.0, 50.0, 200.0]) * R_g,
        "label": r"$z_0/R_g$",
        "filename": "SEDs_z0s.pdf",
        "format": ".3g",
    },
    "gammaIndex": {
        "values": [0.1, 0.3, 0.5],
        "label": r"$\mathrm{gammaIndex}$",
        "filename": "SEDs_gammaIndex.pdf",
        "format": ".2g",
    },
    "geomIndex_out": {
        "values": [0.25, 0.5, 0.75, 1.0],
        "label": r"$\mathrm{geomIndex}_{\rm out}$",
        "filename": "SEDs_geomIndex_out.pdf",
        "format": ".2g",
    },
    "fnT": {
        "values": [0.01, 0.05, 0.1, 0.3],
        "label": r"$f_{\rm nT}$ (with coupled $L_{\rm j}$)",
        "filename": "SEDs_fnT.pdf",
        "format": ".2g",
    },
    "facc": {
        "values": [1e-5, 1e-4, 1e-3, 1e-2],
        "label": r"$\eta_{\rm acc}$",
        "filename": "SEDs_facc.pdf",
        "format": ".1e",
    },
    "p": {
        "values": [2.01, 2.2, 2.4],
        "label": r"$p$",
        "filename": "SEDs_p.pdf",
        "format": ".3g",
    },
    "zmax_inj": {
        "values": np.array([1e-2, 1e-1, 1.0]) * REFERENCE["zmax"],
        "label": r"$z_{\rm max,e}/z_{\rm max}$",
        "filename": "SEDs_zmax_inj.pdf",
        "format": ".2g",
    },
}


def load_jetcode():
    try:
        import JetCode
    except ImportError as error:
        raise SystemExit(
            "JetCode is unavailable. Build it first with 'make -C Code f2py'."
        ) from error
    return JetCode


def configure_model(jetcode, output_dir):
    """Set mutable module state needed by both f2py entry points."""
    sp = jetcode.system_parameters
    sp.set_outdir(str(output_dir))
    sp.z0_rg = REFERENCE["z0"] / R_g
    sp.zmax_inj_rg = REFERENCE["zmax_inj"] / R_g
    sp.gammaindex = REFERENCE["gammaIndex"]
    sp.geomindex_out = REFERENCE["geomIndex_out"]
    sp.fnt = REFERENCE["fnT"]
    sp.facc = REFERENCE["facc"]
    sp.p = REFERENCE["p"]
    sp.derive_parameters()


def run_syn(jetcode, parameters, output_dir):
    configure_model(jetcode, output_dir)
    values = jetcode.syn_model(
        parameters["sigma0"],
        parameters["zmax"],
        parameters["Lj"],
        parameters["a_inj"],
        parameters["zmax_inj"],
        parameters["gammaIndex"],
        parameters["geomIndex_out"],
        parameters["z0"],
        parameters["fnT"],
        parameters["facc"],
        parameters["p"],
        E,
    )
    return np.asarray(values)[: len(E)]


def run_power(jetcode, parameters, output_dir):
    configure_model(jetcode, output_dir)
    result = jetcode.power_model(
        parameters["sigma0"],
        parameters["zmax"],
        parameters["Lj"],
        parameters["a_inj"],
        Z,
    )
    return tuple(np.asarray(component)[: len(Z)] for component in result)


def format_value(parameter, value, fmt):
    display_value = value / R_g if parameter in {"z0", "zmax", "zmax_inj"} else value
    return format(display_value, fmt)


def plot_sed_sweep(parameter, values, seds, output_dir):
    settings = SWEEPS[parameter]
    count = len(values)
    cmap = matplotlib.colormaps["copper_r"].resampled(count)
    plot_colors = cmap(np.linspace(0.0, 1.0, count))
    bounds = np.arange(0.5, count + 1.5)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, (ax, ax_radio) = plt.subplots(1, 2, figsize=(9, 4))
    for index, sed in enumerate(seds):
        ax.loglog(E / h, E * sed * LUMINOSITY_FACTOR, color=plot_colors[index], lw=2)
        ax_radio.loglog(E / h / 1e9, sed * h / 1e-23, color=plot_colors[index], lw=2)

    ax.set_xlim(SED_XMIN, SED_XMAX)
    ax.set_ylim(SED_YMIN, SED_YMAX)
    ax_radio.set_xlim(RADIO_XMIN, RADIO_XMAX)
    ax_radio.set_ylim(RADIO_YMIN, RADIO_YMAX)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax_radio.set_xscale("log")
    ax_radio.set_yscale("log")
    ax.set_xlabel(r"$\nu$ [Hz]")
    ax.set_ylabel(r"$\nu L_\nu$ [erg s$^{-1}$]")
    ax_radio.set_xlabel(r"$\nu$ [GHz]")
    ax_radio.set_ylabel(r"$S_\nu$ [Jy]")
    for axis in (ax, ax_radio):
        axis.grid(which="major", color="0.8")
        axis.grid(which="minor", axis="y", color="0.9", linestyle="--")
        axis.tick_params(direction="in", which="both", top=True, right=True)

    top = ax.twiny()
    top.set_xscale("log")
    top.set_xlim(E[0] / eV, E[-1] / eV)
    top.set_xlabel(r"$E$ [eV]")
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    cbar = fig.colorbar(scalar, ax=ax_radio, ticks=np.arange(1, count + 1), pad=0.02)
    cbar.ax.set_yticklabels(
        [format_value(parameter, value, settings["format"]) for value in values]
    )
    cbar.set_label(settings["label"])
    fig.tight_layout()
    fig.savefig(output_dir / settings["filename"], bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)


def plot_power_sweep(values, components, output_dir):
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = matplotlib.colormaps["copper_r"].resampled(len(values))
    for index, value in enumerate(values):
        color = cmap(index / max(1, len(values) - 1))
        ax.loglog(Z / R_g, components[index][2], color=color, label=f"$a_{{inj}}={value:g}$")
        ax.loglog(Z / R_g, components[index][5], color=color, linestyle=":")
    ax.set_xlabel(r"$z/R_g$")
    ax.set_ylabel(r"Power [erg s$^{-1}$]")
    ax.set_title("Jet power components")
    ax.grid(which="both", color="0.85")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "jet_power_components.pdf", bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)


def run_sweep(name, jetcode, output_dir, scratch_dir):
    settings = SWEEPS[name]
    values = np.asarray(settings["values"], dtype=float)
    results = []
    for index, value in enumerate(values):
        parameters = REFERENCE.copy()
        parameters[name] = value
        if name == "fnT":
            parameters["Lj"] = REFERENCE["Lj"] * REFERENCE["fnT"] / value
        print(f"  {name}={value:g} ({index + 1}/{len(values)})", flush=True)
        results.append(run_syn(jetcode, parameters, scratch_dir))
    plot_sed_sweep(name, values, results, output_dir)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sweep",
        choices=["all", "power", *SWEEPS],
        default="all",
        help="Run one sweep, the power plot, or all experiments (default: all).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "parameter_exploration",
        help="Directory for generated plots (default: diagnostics/parameter_exploration).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    jetcode = load_jetcode()
    with tempfile.TemporaryDirectory(prefix="extendedjet_explore_") as scratch:
        scratch_dir = Path(scratch)
        names = list(SWEEPS) if args.sweep == "all" else ([args.sweep] if args.sweep != "power" else [])
        for name in names:
            print(f"Running {name} sweep")
            run_sweep(name, jetcode, args.output_dir, scratch_dir)

        if args.sweep in {"all", "power"}:
            values = np.asarray([-0.99, 0.01, 1.0])
            components = []
            for value in values:
                parameters = REFERENCE.copy()
                parameters["a_inj"] = value
                print(f"  power a_inj={value:g}", flush=True)
                components.append(run_power(jetcode, parameters, scratch_dir))
            plot_power_sweep(values, components, args.output_dir)

    print(f"Wrote results to {args.output_dir}")


if __name__ == "__main__":
    main()
