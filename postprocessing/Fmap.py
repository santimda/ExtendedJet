"""Build a projected synchrotron luminosity map from a model output.

The map is sampled in the plane of the sky and assumes that each jet segment
emits uniformly across its projected cross-section. Beam convolution is not
part of this script yet; the result is an intrinsic, numerically sampled map.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from postprocessing import plot_utils  # Applies the shared plotting style.
except ImportError:  # Direct execution from the postprocessing directory.
    import plot_utils


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NML = ROOT / "configs" / "m87_reference.nml"
DEFAULT_OUTPUT = ROOT / "output"

# CGS conversion constants.
h = 6.6260755e-27
eV = 1.602176e-12
kpc = 3.0856e21
Mpc = 1.0e3 * kpc


def _read_namelist_values(filename):
    """Read the scalar values needed to project the map from a namelist."""
    text = Path(filename).read_text()
    values = {}
    for name in ("dist_mpc", "incl_deg"):
        match = re.search(
            rf"^\s*{name}\s*=\s*([^,\s/]+)",
            text,
            flags=re.IGNORECASE | re.MULTILINE,
        )
        if match is None:
            raise ValueError(f"{name} was not found in namelist {filename}")
        try:
            values[name] = float(match.group(1).replace("d", "e").replace("D", "E"))
        except ValueError as error:
            raise ValueError(f"Invalid {name} in namelist {filename}") from error
    return values


def read_SEDsyn_cuts(filename):
    """Read the sorted synchrotron SED and infer its grids.

    The file contains photon energy, E*L(E), and segment number. Synchrotron
    emission is written for segments 2 through Ns, so the returned segment
    numbers are also used to select matching rows from quant_FF.dat.
    """
    data = np.loadtxt(filename, comments="#", usecols=(0, 1, 2))
    if data.ndim != 2 or data.shape[1] != 3:
        raise ValueError(f"Expected three columns in {filename}")

    segment_numbers = data[:, 2].astype(int)
    if not np.allclose(segment_numbers, data[:, 2]):
        raise ValueError(f"Non-integer segment number found in {filename}")
    segments = np.unique(segment_numbers)
    if len(segments) == 0:
        raise ValueError(f"No SED data found in {filename}")

    first = data[segment_numbers == segments[0]]
    Esy = first[np.argsort(first[:, 0]), 0]
    Lsy = np.empty((len(Esy), len(segments)))
    for index, segment in enumerate(segments):
        rows = data[segment_numbers == segment]
        rows = rows[np.argsort(rows[:, 0])]
        if len(rows) != len(Esy) or not np.allclose(rows[:, 0], Esy):
            raise ValueError(f"Inconsistent energy grid for segment {segment} in {filename}")
        Lsy[:, index] = rows[:, 1] / Esy

    return Esy, Lsy, segments


def get_Lsyn_z(Esy, Lsy, nu_target=43.0):
    """Interpolate the segment luminosity at ``nu_target`` in GHz."""
    E_nu_target = h * nu_target * 1.0e9 / eV
    if E_nu_target < Esy[0] or E_nu_target > Esy[-1]:
        raise ValueError(
            f"Requested frequency {nu_target:g} GHz is outside the SED range "
            f"{Esy[0] * eV / h / 1.0e9:g}-{Esy[-1] * eV / h / 1.0e9:g} GHz"
        )
    return np.array([np.interp(E_nu_target, Esy, Lsy[:, index]) for index in range(Lsy.shape[1])])


def get_z_rz(filename):
    """Read segment positions and outer radii from ``quant_FF.dat``."""
    data = np.loadtxt(filename, usecols=(0, 2))
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError(f"Expected z and radius columns in {filename}")
    return data[:, 0], data[:, 1]


def make_map(L_z, z, r_z, incl, convert, Nx=512):
    """Project segment luminosities onto a uniformly sampled sky map."""
    if len(L_z) != len(z) or len(z) != len(r_z):
        raise ValueError("SED and jet-geometry segment counts do not match")
    if Nx < 2:
        raise ValueError("Nx must be at least 2")

    y = z * np.sin(incl) * convert
    x = np.linspace(0.0, r_z[-1], Nx) * convert
    dx = x[1] - x[0]
    r_z_proj = r_z * convert
    L_map = np.zeros((len(x), len(y)))

    for j in range(len(y)):
        if r_z_proj[j] <= 0.0:
            continue
        inside = x < r_z_proj[j]
        L_map[inside, j] = (
            2.0 * L_z[j] * np.sqrt(r_z_proj[j] ** 2 - x[inside] ** 2)
            / (np.pi * r_z_proj[j] ** 2)
            * dx
        )
        total = np.sum(L_map[:, j])
        if total > 0.0:
            L_map[:, j] *= 0.5 * L_z[j] / total

    x_map = np.concatenate((-x[::-1], x))
    F_map = np.concatenate((L_map[::-1, :], L_map), axis=0)
    return x_map, y, F_map


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frequency", type=float, default=43.0, help="Target frequency in GHz (default: 43)")
    parser.add_argument("--nml", type=Path, default=DEFAULT_NML, help="Model namelist (default: configs/m87_reference.nml)")
    args = parser.parse_args()

    values = _read_namelist_values(args.nml)
    distance = values["dist_mpc"] * Mpc
    incl = np.deg2rad(values["incl_deg"])
    convert = 206265.0e3 / distance
    flux = 1.0e26 * (h / eV) / (4.0 * np.pi * distance ** 2)

    sed_file = DEFAULT_OUTPUT / "SEDe_syn_cuts_ord.dat"
    quant_file = DEFAULT_OUTPUT / "quant_FF.dat"
    Esy, Lsy, segments = read_SEDsyn_cuts(sed_file)
    L_z = get_Lsyn_z(Esy, Lsy, args.frequency)
    z_all, r_all = get_z_rz(quant_file)

    rows = segments - 1
    if np.any(rows < 0) or np.any(rows >= len(z_all)):
        raise ValueError("SED segment numbers do not match quant_FF.dat")
    z = z_all[rows]
    r_z = r_all[rows]
    x_map, y, F_map = make_map(L_z, z, r_z, incl, convert)
    F_map *= flux

    output = DEFAULT_OUTPUT / f"L_map_{args.frequency:g}GHz.pdf"
    fig, ax = plt.subplots()
    ax.set_xlabel("$x$ [mas]")
    ax.set_ylabel("$y$ [mas]")
    contour = ax.contour(
        y,
        x_map,
        np.log10(np.where(F_map > 0.0, F_map, np.nan)),
        levels=256,
        cmap="afmhot_r",
    )
    fig.colorbar(contour, ax=ax, label="log10 flux [mJy]")
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
