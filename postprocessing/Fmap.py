"""Create two beam-convolved jet intensity maps and one spectral-index map.

The default invocation creates 22 and 43 GHz maps with 1.2 and 0.31 mas
circular beams.  The spectral-index map uses both frequencies convolved to the
lower-frequency beam, as required for a resolved spectral-index measurement.

Examples::

    python postprocessing/Fmap.py
    python postprocessing/Fmap.py --frequencies 22 43 --beams 1.2 0.31
    python postprocessing/Fmap.py --beams 1.2,1.0 0.31,0.25

Beam values are FWHM in mas.  A scalar is circular; ``major,minor`` is
elliptical.  Beam position angles are intentionally not supported yet.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import fftconvolve

try:
    from postprocessing import plot_utils  # noqa: F401
except ImportError:
    import plot_utils  # noqa: F401


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NML = ROOT / "configs" / "m87_reference.nml"
MAS = np.deg2rad(1.0 / 3.6e6)
MPC = 3.0856775814913673e24
H = 6.6260755e-27
GHZ = 1.0e9
PIXEL_MAS = 0.004
HALF_WIDTH_MAS = 2.5
DEFAULT_FREQUENCIES = (22.0, 43.0)
DEFAULT_BEAMS = ((1.2, 1.2), (0.31, 0.31))
DEFAULT_Y_MAX = {22.0: 8.0, 43.0: 4.0}
JET_PA_DEG = -72.0


def _value(text, name):
    match = re.search(rf"^\s*{name}\s*=\s*([^,\n]+)", text, re.IGNORECASE | re.MULTILINE)
    if match is None:
        raise ValueError(f"{name} was not found in the namelist")
    return match.group(1).strip().strip("'\"")


def _read_namelist(filename):
    text = Path(filename).read_text()
    values = {}
    for name in ("dist_mpc", "incl_deg", "M_BH"):
        raw = _value(text, name)
        values[name] = float(raw.replace("d", "e").replace("D", "E"))
    values["outdir"] = _value(text, "outdir")
    raw = _value(text, "grav_redshift")
    values["grav_redshift"] = raw.lower() in {".true.", "true", "on"}
    return values


def _resolve_outdir(nml, value):
    path = Path(value)
    if path.is_absolute():
        return path
    # The documented Fortran command runs from Code/, so relative outdir
    # values in a namelist have that same meaning for post-processing.
    return (ROOT / "Code" / path).resolve()


def _parse_beam(value):
    parts = value.strip().strip("()").replace("x", ",").split(",")
    try:
        numbers = tuple(float(part) for part in parts if part.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"invalid beam {value!r}") from error
    if len(numbers) == 1:
        numbers = (numbers[0], numbers[0])
    if len(numbers) != 2 or min(numbers) <= 0:
        raise argparse.ArgumentTypeError("a beam must be one positive value or major,minor")
    return numbers


def _load_segments(path, electron_path, grav_redshift, m_bh):
    quantities = np.loadtxt(path)
    if quantities.ndim != 2 or quantities.shape[1] < 14:
        raise ValueError(f"Expected at least 14 columns in {path}")
    electron = np.loadtxt(electron_path)
    if electron.ndim != 2 or electron.shape[1] < 3:
        raise ValueError(f"Expected three columns in {electron_path}")

    ns = len(quantities)
    segment_ids = electron[:, 2].astype(int)
    if not np.array_equal(np.unique(segment_ids), np.arange(1, ns + 1)):
        raise ValueError(f"Electron segments do not match {path}")
    first = electron[segment_ids == 1]
    energies = first[np.argsort(first[:, 0]), 0]
    spectra = np.empty((len(energies), ns))
    for index in range(ns):
        rows = electron[segment_ids == index + 1]
        rows = rows[np.argsort(rows[:, 0])]
        if len(rows) != len(energies) or not np.allclose(rows[:, 0], energies):
            raise ValueError(f"Inconsistent electron energy grid in {electron_path}")
        spectra[:, index] = rows[:, 1]

    z = np.abs(quantities[:, 0])
    radii = np.sqrt(z * z + 0.5 * (quantities[:, 2] ** 2 + quantities[:, 3] ** 2))
    rg = 6.67430e-8 * (m_bh * 1.98847e33) / (2.99792458e10 ** 2)
    if grav_redshift:
        grav = np.sqrt(np.clip(1.0 - 2.0 * rg / radii, 0.0, None))
    else:
        grav = np.ones(ns)
    return {
        "z": z,
        "B": quantities[:, 1],
        "R": quantities[:, 2],
        "R_in": quantities[:, 3],
        "delta": quantities[:, 5],
        "dz": quantities[:, 10],
        "g": grav,
        "E": energies,
        "N": spectra,
    }


def _load_jetcode(nml):
    sys.path.insert(0, str(ROOT / "Code"))
    try:
        import JetCode
    except ImportError as error:
        raise RuntimeError("Build the f2py module first with `make -C Code f2py`") from error
    if not hasattr(JetCode, "local_synchrotron"):
        raise RuntimeError("JetCode is stale; rebuild it with `make -C Code f2py`")
    JetCode.system_parameters.read_namelist(str(nml))
    JetCode.system_parameters.derive_parameters()
    return JetCode


def _emission(JetCode, segments, frequency):
    j_nu = np.zeros(len(segments["z"]))
    alpha = np.zeros_like(j_nu)
    for index in range(len(j_nu)):
        if not np.any(segments["N"][:, index] > 0):
            continue
        doppler = segments["delta"][index] * segments["g"][index]
        if doppler <= 0:
            continue
        j_nu[index], alpha[index] = JetCode.local_synchrotron(
            H * frequency * GHZ / doppler,
            segments["B"][index],
            segments["N"][:, index],
            segments["E"],
        )
    return j_nu, alpha


def _jet_tables(JetCode, outdir, nml_values, frequencies):
    tables = {}
    electron_path = outdir / "Ne_cuts.dat"
    for frequency in frequencies:
        tables[frequency] = {}
        for suffix, side in (("", 1), ("CJ", -1)):
            segments = _load_segments(
                outdir / f"quant_FF{suffix}.dat",
                electron_path,
                nml_values["grav_redshift"],
                nml_values["M_BH"],
            )
            j_nu, alpha = _emission(JetCode, segments, frequency)
            good = (j_nu > 0) & (alpha >= 0) & (segments["z"] > 0)
            if not np.any(good):
                raise ValueError(f"No emitting segments found at {frequency:g} GHz")
            tables[frequency][side] = {
                "log_z": np.log(segments["z"][good]),
                "log_j": np.log(j_nu[good]),
                "log_alpha": np.log(np.maximum(alpha[good], np.finfo(float).tiny)),
                "doppler": segments["delta"][good] * segments["g"][good],
                "R": segments["R"][good],
                "R_in": segments["R_in"][good],
            }
    return tables


def _chord(x, outer, inner):
    outer_chord = 2.0 * np.sqrt(np.clip(outer[:, None] ** 2 - x[None, :] ** 2, 0.0, None))
    inner_chord = 2.0 * np.sqrt(np.clip(inner[:, None] ** 2 - x[None, :] ** 2, 0.0, None))
    return outer_chord - inner_chord


def jet_intensity(tables, frequency, x_mas, y_mas, distance, inclination):
    image = np.zeros((len(y_mas), len(x_mas)))
    x = x_mas * MAS * distance
    for side in (1, -1):
        table = tables[frequency][side]
        rows = np.where(side * y_mas > 0)[0]
        z = np.abs(y_mas[rows]) * MAS * distance / np.sin(inclination)
        inside = (np.log(z) >= table["log_z"][0]) & (np.log(z) <= table["log_z"][-1])
        rows, z = rows[inside], z[inside]
        if len(rows) == 0:
            continue
        interp = lambda values: np.interp(np.log(z), table["log_z"], values)
        j = np.exp(interp(table["log_j"]))
        alpha = np.exp(interp(table["log_alpha"]))
        doppler = np.interp(np.log(z), table["log_z"], table["doppler"])
        radius = np.exp(interp(np.log(table["R"])))
        inner = np.exp(interp(np.log(np.maximum(table["R_in"], np.finfo(float).tiny))) )
        chord = _chord(x, radius, inner)
        tau = alpha[:, None] * chord / (doppler[:, None] * np.sin(inclination))
        image[rows] = (doppler ** 3 * j / alpha)[:, None] * (-np.expm1(-tau))
    return image


def convolve(image, beam, pixel_mas):
    sy, sx = (value / (2.0 * np.sqrt(2.0 * np.log(2.0))) for value in beam)
    ky = np.arange(-int(4 * sy / pixel_mas), int(4 * sy / pixel_mas) + 1) * pixel_mas
    kx = np.arange(-int(4 * sx / pixel_mas), int(4 * sx / pixel_mas) + 1) * pixel_mas
    kernel = np.exp(-0.5 * (ky[:, None] ** 2 / sy ** 2 + kx[None, :] ** 2 / sx ** 2))
    pixel_flux = image * (pixel_mas * MAS) ** 2 * 1e26
    return fftconvolve(pixel_flux, kernel, mode="same")


def _sky(x, y, jet_pa):
    pa = np.deg2rad(jet_pa)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    ra = yy * np.sin(pa) + xx * np.cos(pa)
    dec = yy * np.cos(pa) - xx * np.sin(pa)
    return -ra, dec


def _beam_patch(ax, beam, jet_pa):
    """Draw the beam in the lower-left sky corner.

    Beam major/minor axes are defined along/across the jet, so the major-axis
    position angle is the plotted jet position angle until beam PA is exposed
    as a separate command-line parameter.
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    margin = 0.08
    x = xlim[0] + margin * (xlim[1] - xlim[0]) + beam[1] / 2.0
    y = ylim[0] + margin * (ylim[1] - ylim[0]) + beam[0] / 2.0
    # Matplotlib measures angles from the mirrored RA axis; jet_pa is measured
    # east of north, hence the 90-degree conversion plus the reflection.
    ax.add_patch(Ellipse((x, y), beam[1], beam[0], angle=90.0 + jet_pa, facecolor="grey", edgecolor="none", zorder=10))


def _add_colorbar(fig, ax, plot, label):
    divider = make_axes_locatable(ax)
    colorbar_axis = divider.append_axes("right", size="4%", pad=0.08)
    fig.colorbar(plot, cax=colorbar_axis, label=label)


def _plot_intensity(x, y, image, frequency, beam, jet_pa, output):
    ra, dec = _sky(x, y, jet_pa)
    maximum = float(np.max(image))
    if not np.isfinite(maximum) or maximum <= 0.0:
        raise ValueError(f"The {frequency:g} GHz map has no positive intensity")
    lower = maximum * 1e-5
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    plot = ax.pcolormesh(ra, dec, np.maximum(image, lower), norm=LogNorm(lower, maximum), cmap="afmhot_r", shading="auto", rasterized=True)
    levels = []
    level = 0.9 * maximum
    while level >= lower:
        levels.append(level)
        level /= 3.0
    ax.contour(ra, dec, image, levels=levels[::-1], colors="grey", linewidths=0.45)
    _add_colorbar(fig, ax, plot, rf"$S_{{{frequency:g}}}$ [mJy beam$^{{-1}}$]")
    ax.set_xlabel("relative RA [mas]")
    ax.set_ylabel("relative Dec [mas]")
    ax.set_aspect("equal")
    _beam_patch(ax, beam, jet_pa)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_spectral_index(x, y, low_image, high_image, low, high, beam, jet_pa, output):
    with np.errstate(divide="ignore", invalid="ignore"):
        spectral_index = np.log(high_image / low_image) / np.log(high / low)
    low_threshold = float(np.max(low_image)) * 1e-5
    high_threshold = float(np.max(high_image)) * 1e-5
    spectral_index = np.where((low_image >= low_threshold) & (high_image >= high_threshold), spectral_index, np.nan)
    ra, dec = _sky(x, y, jet_pa)
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    plot = ax.pcolormesh(ra, dec, spectral_index, cmap="jet", norm=Normalize(-1.5, 0.5), shading="auto", rasterized=True)
    _add_colorbar(fig, ax, plot, r"spectral index $\alpha$")
    ax.set_title(rf"{low:g}--{high:g} GHz")
    ax.set_xlabel("relative RA [mas]")
    ax.set_ylabel("relative Dec [mas]")
    ax.set_aspect("equal")
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nml", type=Path, default=DEFAULT_NML)
    parser.add_argument("--frequencies", nargs=2, type=float, default=DEFAULT_FREQUENCIES, metavar=("LOW", "HIGH"))
    parser.add_argument("--beams", nargs=2, type=_parse_beam, default=DEFAULT_BEAMS, metavar=("LOW_BEAM", "HIGH_BEAM"), help="one scalar or major,minor FWHM in mas per frequency")
    parser.add_argument("--jet-pa", type=float, default=JET_PA_DEG, help="approaching-jet position angle in degrees")
    parser.add_argument("--outdir", type=Path, help="override the namelist output directory")
    args = parser.parse_args()
    low, high = args.frequencies
    if not 0 < low < high:
        parser.error("frequencies must be positive and strictly increasing")
    nml_values = _read_namelist(args.nml)
    outdir = args.outdir.resolve() if args.outdir else _resolve_outdir(args.nml, nml_values["outdir"])
    figure_dir = outdir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    distance = nml_values["dist_mpc"] * MPC
    inclination = np.deg2rad(nml_values["incl_deg"])
    JetCode = _load_jetcode(args.nml)
    tables = _jet_tables(JetCode, outdir, nml_values, (low, high))
    x = np.arange(-HALF_WIDTH_MAS, HALF_WIDTH_MAS + PIXEL_MAS / 2, PIXEL_MAS)
    ymax = max(DEFAULT_Y_MAX.get(low, 8.0), DEFAULT_Y_MAX.get(high, 4.0))
    y = np.arange(-ymax, ymax + PIXEL_MAS / 2, PIXEL_MAS)
    intrinsic = {frequency: jet_intensity(tables, frequency, x, y, distance, inclination) for frequency in (low, high)}
    low_map = convolve(intrinsic[low], args.beams[0], PIXEL_MAS)
    high_map = convolve(intrinsic[high], args.beams[1], PIXEL_MAS)
    high_low_beam = convolve(intrinsic[high], args.beams[0], PIXEL_MAS)
    low_output = figure_dir / f"Map_I_{low:g}GHz.pdf"
    high_output = figure_dir / f"Map_I_{high:g}GHz.pdf"
    spectral_output = figure_dir / "Map_spix.pdf"
    _plot_intensity(x, y, low_map, low, args.beams[0], args.jet_pa, low_output)
    _plot_intensity(x, y, high_map, high, args.beams[1], args.jet_pa, high_output)
    _plot_spectral_index(x, y, low_map, high_low_beam, low, high, args.beams[0], args.jet_pa, spectral_output)
    print(json.dumps({"outdir": str(figure_dir), "frequencies_GHz": [low, high], "beams_mas": args.beams, "jet_pa_deg": args.jet_pa}, indent=2))


if __name__ == "__main__":
    main()
