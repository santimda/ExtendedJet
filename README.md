# ExtendedJet

Non-thermal multi-zone emission model for relativistic jets. The model is partly based on Molina, del Palacio & Bosch Ramon 2019 (ADS link: https://ui.adsabs.harvard.edu/abs/2019A%26A...629A.129M/abstract). The code is written as Fortran 90 modules, compiled with Make, and run as `./run.x params.nml`. Build products can be safely removed using `make clean` from `Code/`.

Some further details of each code are given below.

WARNING: The shell variable LC_NUMERIC="en_US.UTF-8" must be set. The scripts ".sh" are used to reorder some files and it is important that the bash decimal separator is a "." and not a ",". Do 'export LC_NUMERIC="en_US.UTF-8"' before running or the SED might be unsorted.

## Environment

Create the minimal environment for building and running the model, f2py, tests,
and standard post-processing with mamba:

```bash
mamba env create -f environment.yml
mamba activate extendedjet
```

For the optional notebook dependencies, use the extended environment instead:

```bash
mamba env create -f environment-notebooks.yml
mamba activate extendedjet-notebooks
```

The Python plotting scripts use LaTeX by default through Matplotlib's
`text.usetex`. Install the system package `texlive-latex-base` to use that
default. If LaTeX is unavailable, set `text.usetex` to `False` in
`postprocessing/plot_utils.py`. The extended environment adds Jupyter, SciPy,
Astropy, pandas, `agnpy`, `naima`, Gnuplot, and `smplotlib`. `agnpy` is used by
`emisividad1.ipynb` for synchrotron/absorption comparisons; `naima` is used by
`Syn_SSC_sed.ipynb` for synchrotron and SSC comparisons. Neither is required
by the Fortran model or its f2py wrapper. External RIAF data used by some
legacy notebooks are not included here.

## Layout

- `Code/`: Fortran sources, Makefile, and generated build files.
- `configs/`: runtime namelists.
- `data/`: legacy and observational input data.
- `diagnostics/`: diagnostic and validation notebooks and plotting scripts.
- `figures/`: repository figures and diagrams.
- `postprocessing/`: synthetic-map notebooks, Python utilities, and data.
- `tests/`: regression and f2py checks.
- `validation/`: optional physics-validation notebooks.
- `output/`: generated reference-model outputs.
- `tests/output/`: disposable outputs from the reduced test configuration.

Run the model from the repository root with a configuration file:

```text
make -C Code
make -C Code test-unit
(cd Code && make test)
(cd Code && OMP_NUM_THREADS=1 ./run.x ../configs/test.nml)
(cd Code && OMP_NUM_THREADS=1 ./run.x ../configs/m87_reference.nml)
```

The reduced run writes to `tests/output/`; the reference run writes to
`output/`. Both directories are disposable and ignored by Git. Diagnostic
Gnuplot scripts are run from `diagnostics/`, and post-processing utilities are
run from `postprocessing/`.

The reduced regression configuration uses `Ns=100` and `numbins=250`, a deliberately lighter but still resolved grid
for fast test execution. The production M87 configuration uses the converged resolution described below.

The notebook parameter exploration is also available as a script. Build the f2py module and run all configured
sweeps with:

```bash
make -C Code f2py
OMP_NUM_THREADS=1 python diagnostics/explore_parameters.py
```

Plots are written to `diagnostics/parameter_exploration/` by default. Use `--sweep NAME` to run one SED sweep or
`--sweep power` for the injection power-components plot; `--output-dir PATH` changes the destination. Intermediate
Fortran files are written to a temporary directory and removed automatically. The script uses the notebook reference
dictionary and runs calls sequentially because the Fortran modules retain mutable state.

`make -C Code test` runs the numerical unit tests and a reduced Fortran
regression. The regression runs `run.x` with `configs/test.nml` in a temporary
output directory and compares the jet and counter-jet fluxes with the frozen
references in `tests/reference/test/` over the radio band.

The current model uses only the cosmic microwave background (CMB) as its
external photon field. The generic IC kernel remains capable of receiving
other black-body-like fields, but no stellar or other external photon field is
configured in this version.

The production M87 grid uses `Ns=700` jet segments and `numbins=1100` particle-transport energy bins. These were
selected from finer full-model convergence tests using a mixed tolerance of 3% for each aperture flux and 1% for total
flux. The compile-time electron output grid remains `nIn=150`; increasing it to 300 changed the fluxes by less than
1%. High-resolution reference runs used an isolated `numdt=300000` build. The selected `Ns=700` grid remains within
the production `numdt=100000` transport-history limit; still finer runs can exceed it and should not be used as
convergence references without also revisiting that limit.

## Important files

**Makefile**: production, f2py, and unit-test build orchestration.

**Makefile.depend**: details which modules depend on which other modules (needed for compilation)

**constants.f90**: useful physical and numerical constants. It is imported by **global.f90** to be used by all other modules as well.

**numerical_utils.f90**: reusable grids, integration, interpolation, safe logarithms, and linear-algebra utilities.

**system_parameters.f90**: parameters of the specific system being studied. It is imported by **globals.f90** to be used by all other modules as well.

**globals.f90**: simulation grid sizes, shared state arrays, and common radiation functions.

**controls.f90**: compile-time and namelist-controlled calculation switches.

**trajectory.f90**: calculates the jet or counter-jet trajectory (trivial for a vertical jet).

**evolution.f90**: calculates the evolution of thermodynamical quantities along the jet.

**losses.f90**: contains functions for calculating energy losses for electrons.

**edistribution.f90**: calculates the non-thermal electron population along the jet and saves it in a plain text file.

!**XXrad.f90**: calculates the SED by the respective process reading the corresponding non-thermal particle population from the file generated by **edistribution.f90**.
