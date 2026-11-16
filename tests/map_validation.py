"""Validate the reference M87 map flux against the Fortran calculation."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "postprocessing"))
import Fmap  # noqa: E402


def fortran_flux_mjy(JetCode, outdir, nml_values, frequency, distance, inclination, ymax):
    total = 0.0
    for suffix in ("", "CJ"):
        segments = Fmap._load_segments(
            outdir / f"quant_FF{suffix}.dat",
            outdir / "Ne_cuts.dat",
            nml_values["grav_redshift"],
            nml_values["M_BH"],
        )
        projected_z = segments["z"] * np.sin(inclination) / (Fmap.MAS * distance)
        for index in np.where(projected_z <= ymax)[0]:
            doppler = segments["delta"][index] * segments["g"][index]
            if doppler <= 0.0:
                continue
            flux_per_energy = JetCode.local_dflux(
                segments["N"][:, index],
                segments["B"][index],
                Fmap.H * frequency * Fmap.GHZ / doppler,
                segments["E"],
                segments["delta"][index],
                segments["g"][index],
                segments["R"][index],
                segments["R_in"][index],
            )
            total += flux_per_energy * Fmap.H * segments["dz"][index] * 1e26
    return float(total)


def main():
    nml = Fmap.DEFAULT_NML
    nml_values = Fmap._read_namelist(nml)
    outdir = Fmap._resolve_outdir(nml, nml_values["outdir"])
    required = [outdir / "Ne_cuts.dat", outdir / "quant_FF.dat", outdir / "quant_FFCJ.dat"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Reference map output is missing; run configs/m87_reference.nml first: " + ", ".join(missing))

    distance = nml_values["dist_mpc"] * Fmap.MPC
    inclination = np.deg2rad(nml_values["incl_deg"])
    frequencies = (22.0, 43.0)
    JetCode = Fmap._load_jetcode(nml)
    tables = Fmap._jet_tables(JetCode, outdir, nml_values, frequencies)
    x = np.arange(-Fmap.HALF_WIDTH_MAS, Fmap.HALF_WIDTH_MAS + Fmap.PIXEL_MAS / 2, Fmap.PIXEL_MAS)
    ymax = max(Fmap.DEFAULT_Y_MAX.values())
    y = np.arange(-ymax, ymax + Fmap.PIXEL_MAS / 2, Fmap.PIXEL_MAS)

    result = {}
    for frequency in frequencies:
        image = Fmap.jet_intensity(tables, frequency, x, y, distance, inclination)
        map_flux = float(np.sum(image) * (Fmap.PIXEL_MAS * Fmap.MAS) ** 2 * 1e26)
        code_flux = fortran_flux_mjy(JetCode, outdir, nml_values, frequency, distance, inclination, ymax)
        relative_difference = abs(map_flux - code_flux) / max(abs(code_flux), np.finfo(float).tiny)
        result[f"{frequency:g}GHz"] = {
            "map_mJy": map_flux,
            "fortran_mJy": code_flux,
            "relative_difference": relative_difference,
        }
        if relative_difference > 0.01:
            raise SystemExit(f"{frequency:g} GHz map flux differs from Fortran by {relative_difference:.3%}")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
