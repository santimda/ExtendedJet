"""Generate synthetic aperture fluxes and run the reusable SED fitter."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fitting.config import FluxData, mJy, save_fluxes
from fitting.fit import fit_sed
from fitting.model import ExtendedJetModel


DEFAULT_NML = ROOT / "configs" / "m87_reference.nml"
DEFAULT_OUTDIR = ROOT / "validation" / "output" / "mcmc_validation"
FREQUENCIES_GHZ = np.array([3.0, 10.0, 22.0, 43.0, 86.0])
APERTURES_RG = np.array([30.0, 100.0, 200.0])
FIT_PARAMETERS = ("geomIndex_out", "a_inj")
ERROR_FRACTION = 0.15


def main():
    DEFAULT_OUTDIR.mkdir(parents=True, exist_ok=True)
    model = ExtendedJetModel(DEFAULT_NML, DEFAULT_OUTDIR / "exact_model", ns=40, numbins=100)
    truth = {"Lj": 1e42, "a_inj": 0.01, "geomIndex_out": 0.56}
    fluxes = model(FREQUENCIES_GHZ, APERTURES_RG, truth)
    errors = ERROR_FRACTION * fluxes
    data = FluxData.from_mjy(FREQUENCIES_GHZ, APERTURES_RG, fluxes / mJy, errors / mJy, FIT_PARAMETERS)
    save_fluxes(DEFAULT_OUTDIR / "synthetic_fluxes.json", data)
    result = fit_sed(
        model=model,
        data=data,
        output_dir=DEFAULT_OUTDIR,
        mode="overwrite",
        fixed_parameters=truth,
        truth=truth,
        nwalkers=8,
        nsteps=200,
        seed=1,
    )
    posterior = result.posterior
    for name in FIT_PARAMETERS:
        values = np.asarray(posterior[name], dtype=float)
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"non-finite posterior samples for {name}")
        if not np.quantile(values, 0.05) <= truth[name] <= np.quantile(values, 0.95):
            raise RuntimeError(f"injected value for {name} is outside the 90% posterior interval")


if __name__ == "__main__":
    main()
