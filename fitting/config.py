"""Configuration and flux-data loading for source fits."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

Jy = 1e-23
mJy = 1e-26


@dataclass(frozen=True)
class FluxData:
    frequencies_ghz: np.ndarray
    apertures_rg: np.ndarray
    fluxes_cgs: np.ndarray
    errors_cgs: np.ndarray
    fit_parameters: tuple[str, ...] = ()

    @classmethod
    def from_mjy(cls, frequencies_ghz, apertures_rg, fluxes_mjy, errors_mjy, fit_parameters):
        frequencies = np.asarray(frequencies_ghz, dtype=float)
        apertures = np.asarray(apertures_rg, dtype=float)
        fluxes = _numeric_matrix(fluxes_mjy, "fluxes_mJy")
        errors = _numeric_matrix(errors_mjy, "errors_mJy")
        expected = (len(frequencies), len(apertures))
        if fluxes.shape != expected:
            raise ValueError(f"fluxes_mJy must have shape {expected}")
        if errors is not None and errors.shape != expected:
            raise ValueError(f"errors_mJy must have shape {expected}")
        if frequencies.ndim != 1 or np.any(frequencies <= 0):
            raise ValueError("frequencies_GHz must be one-dimensional and positive")
        if apertures.ndim != 1 or np.any(apertures <= 0):
            raise ValueError("apertures_Rg must be one-dimensional and positive")
        if np.any(np.diff(frequencies) <= 0) or np.any(np.diff(apertures) <= 0):
            raise ValueError("frequencies_GHz and apertures_Rg must be strictly sorted")
        observed = np.isfinite(fluxes)
        if not np.any(observed):
            raise ValueError("at least one flux measurement is required")
        if np.any(fluxes[observed] <= 0):
            raise ValueError("measured fluxes must be positive; use null for missing values")
        if errors is not None:
            if np.any(np.isfinite(errors) != observed):
                raise ValueError("errors_mJy must use null in exactly the missing-flux cells")
            if np.any(errors[observed] <= 0):
                raise ValueError("measured errors must be positive")
        if not fit_parameters:
            raise ValueError("fit_parameters must contain at least one parameter")
        return cls(frequencies, apertures, fluxes * mJy, errors * mJy, tuple(fit_parameters))

    @property
    def observed(self):
        return np.isfinite(self.fluxes_cgs)

    @property
    def fluxes_mjy(self):
        return self.fluxes_cgs / mJy

    @property
    def errors_mjy(self):
        return self.errors_cgs / mJy


def _numeric_matrix(values, name):
    matrix = np.asarray(values, dtype=object)
    result = np.full(matrix.shape, np.nan, dtype=float)
    for index, value in np.ndenumerate(matrix):
        if value is not None:
            try:
                result[index] = float(value)
            except (TypeError, ValueError) as error:
                raise ValueError(f"{name} contains a non-numeric value at {index}") from error
            if not np.isfinite(result[index]):
                raise ValueError(f"{name} contains NaN or infinity at {index}; use null for missing values")
    return result


def load_fluxes(path):
    with Path(path).open() as stream:
        content = json.load(stream)
    required = ("fit_parameters", "frequencies_GHz", "apertures_Rg", "fluxes_mJy", "errors_mJy")
    missing = [name for name in required if name not in content]
    if missing:
        raise ValueError(f"missing flux configuration keys: {missing}")
    return FluxData.from_mjy(
        content["frequencies_GHz"],
        content["apertures_Rg"],
        content["fluxes_mJy"],
        content["errors_mJy"],
        content["fit_parameters"],
    )


def save_fluxes(path, data):
    content = {
        "frequencies_GHz": data.frequencies_ghz.tolist(),
        "apertures_Rg": data.apertures_rg.tolist(),
        "fit_parameters": list(data.fit_parameters),
        "fluxes_mJy": [
            [None if not np.isfinite(value) else float(value) for value in row]
            for row in data.fluxes_mjy
        ],
    }
    content["errors_mJy"] = [
        [None if not np.isfinite(value) else float(value) for value in row]
        for row in data.errors_mjy
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        json.dump(content, stream, indent=2)
        stream.write("\n")
