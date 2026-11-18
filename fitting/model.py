"""Exact-frequency ExtendedJet model interface."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "Code"


class ExtendedJetModel:
    """Evaluate circular-aperture synchrotron flux at requested frequencies."""

    log_parameters = {
        "log_Lj": "Lj",
        "epsilon_e": "fnT",
        "z_0": "z0",
        "z_inj": "zmax_inj",
        "log_sigma_0": "sigma_0",
        "eta_acc": "facc",
    }
    parameter_aliases = {
        "g": "gammaIndex",
    }
    sample_parameter_aliases = {
        "sigma_0": "log_sigma_0",
    }
    parameter_names = {
        "sigma_0", "zmax", "Lj", "a_inj", "zmax_inj", "gammaIndex",
        "geomIndex_out", "z0", "fnT", "facc", "p",
    }

    def __init__(self, namelist, output_dir, ns=None, numbins=None):
        sys.path.insert(0, str(CODE))
        import JetCode

        self.jetcode = JetCode
        self.system_parameters = JetCode.system_parameters
        self.system_parameters.read_namelist(str(Path(namelist).resolve()))
        if ns is not None:
            self.system_parameters.ns = int(ns)
        if numbins is not None:
            self.system_parameters.numbins = int(numbins)
        self.system_parameters.derive_parameters()
        # Fitting apertures are defined as projected circular apertures.
        JetCode.controls.aperture_area = True
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.source_name = Path(namelist).stem
        self.h = float(JetCode.constants.h)
        self.parameters = {
            "sigma_0": float(self.system_parameters.sigma_0),
            "zmax": float(self.system_parameters.zmax),
            "Lj": float(self.system_parameters.lj),
            "a_inj": float(self.system_parameters.a_inj),
            "zmax_inj": float(self.system_parameters.zmax_inj),
            "gammaIndex": float(self.system_parameters.gammaindex),
            "geomIndex_out": float(self.system_parameters.geomindex_out),
            "z0": float(self.system_parameters.z0),
            "fnT": float(self.system_parameters.fnt),
            "facc": float(self.system_parameters.facc),
            "p": float(self.system_parameters.p),
        }

    def __call__(self, frequencies_ghz, apertures_rg, parameters=None, output_dir=None):
        parameters = dict(parameters or {})
        accepted_names = self.parameter_names | set(self.log_parameters) | set(self.parameter_aliases)
        unknown = set(parameters) - accepted_names
        if unknown:
            raise ValueError(f"Unknown model parameters: {sorted(unknown)}")
        for public_name, internal_name in self.log_parameters.items():
            if internal_name in parameters:
                raise ValueError(f"provide either {public_name} or {internal_name}, not both")
        for public_name, internal_name in self.parameter_aliases.items():
            if internal_name in parameters:
                raise ValueError(f"provide either {public_name} or {internal_name}, not both")
        values = self.parameters.copy()
        for name, value in parameters.items():
            if name in self.log_parameters:
                internal_name = self.log_parameters[name]
                value = 10.0 ** float(value)
                if name in {"z_0", "z_inj"}:
                    value *= float(self.system_parameters.r_g)
                values[internal_name] = value
            elif name in self.parameter_aliases:
                values[self.parameter_aliases[name]] = value
            else:
                values[name] = value
        frequencies_ghz = np.asarray(frequencies_ghz, dtype=float)
        apertures_rg = np.asarray(apertures_rg, dtype=float)
        if frequencies_ghz.ndim != 1 or np.any(frequencies_ghz <= 0):
            raise ValueError("frequencies_ghz must be one-dimensional and positive")
        if apertures_rg.ndim != 1 or np.any(apertures_rg <= 0):
            raise ValueError("apertures_rg must be one-dimensional and positive")
        if np.any(np.diff(frequencies_ghz) <= 0) or np.any(np.diff(apertures_rg) <= 0):
            raise ValueError("frequencies_ghz and apertures_rg must be strictly sorted")
        if len(frequencies_ghz) > 300 or len(apertures_rg) > 10:
            raise ValueError("the exact evaluator supports at most 300 frequencies and 10 apertures")
        target_dir = Path(output_dir) if output_dir is not None else self.output_dir
        target_dir.mkdir(parents=True, exist_ok=True)
        self.system_parameters.set_outdir(str(target_dir.resolve()))
        raw = self.jetcode.syn_model_apertures(
            values["sigma_0"], values["zmax"], values["Lj"], values["a_inj"],
            values["zmax_inj"], values["gammaIndex"], values["geomIndex_out"],
            values["z0"], values["fnT"], values["facc"], values["p"],
            self.h * frequencies_ghz * 1e9, apertures_rg,
        )
        result = np.asarray(raw)[:len(frequencies_ghz), :len(apertures_rg)]
        if result.shape != (len(frequencies_ghz), len(apertures_rg)) or not np.all(np.isfinite(result)):
            raise ValueError("exact evaluator returned an invalid flux array")
        return result
