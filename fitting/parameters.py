"""Common fit parameters and universal priors."""

from __future__ import annotations

import bilby


DEFAULT_FIT_PARAMETERS = ("geomIndex_out", "a_inj")
COMMON_PRIORS = {
    "log_Lj": (40.0, 43.0, r"$\log_{10}(L_j/[\mathrm{erg\ s^{-1}}])$"),
    "geomIndex_out": (0.0, 1.0, r"$w$"),
    "a_inj": (-0.5, 2.0, r"$a_e$"),
}


def build_priors(fit_parameters):
    fit_parameters = tuple(fit_parameters)
    if "Lj" in fit_parameters:
        raise ValueError("fit log_Lj instead of Lj")
    unknown = set(fit_parameters) - set(COMMON_PRIORS)
    if unknown:
        raise ValueError(f"No common prior is defined for: {sorted(unknown)}")
    priors = bilby.core.prior.PriorDict()
    for name in fit_parameters:
        minimum, maximum, label = COMMON_PRIORS[name]
        priors[name] = bilby.core.prior.Uniform(minimum, maximum, name=name, latex_label=label)
    return priors
