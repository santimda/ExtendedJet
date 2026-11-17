#!/usr/bin/env python3
"""Fit an ExtendedJet source from a namelist and flux JSON file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fitting.config import load_fluxes
from fitting.parameters import DEFAULT_FIT_PARAMETERS
from fitting.runner import fit_sed
from validation.extendedjet_model import ExtendedJetModel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nml", type=Path, required=True)
    parser.add_argument("--fluxes", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fit-parameters", nargs="+", default=DEFAULT_FIT_PARAMETERS)
    parser.add_argument("--ns", type=int)
    parser.add_argument("--numbins", type=int)
    parser.add_argument("--nwalkers", type=int, default=8)
    parser.add_argument("--nsteps", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    data = load_fluxes(args.fluxes)
    model = ExtendedJetModel(args.nml, args.outdir / "exact_model", ns=args.ns, numbins=args.numbins)
    fit_sed(
        model=model,
        data=data,
        fit_parameters=args.fit_parameters,
        output_dir=args.outdir,
        nwalkers=args.nwalkers,
        nsteps=args.nsteps,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
