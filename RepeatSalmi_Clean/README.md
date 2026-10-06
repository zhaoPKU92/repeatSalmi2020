# Salmi 2020 calculation code

This folder contains only the C++ solver used for the current Figure 1, 3 and 8 calculations, its local headers, build/run/plot scripts, and physical parameters and mass coordinates in `code/model_config.py`. The C++ physics files are unchanged. Results, article files, executables, source archives, tests, diagnostics programs, reference curves and obsolete scripts have been removed from this copy. The original `salmi2020` folder is untouched.

## Build and run

Use Python 3.10 or later and GCC with C++14 and OpenMP. From `code`:

```powershell
py -3.10 -m pip install -r requirements.txt
py -3.10 build.py --target atmosphere
py -3.10 build.py --target figure1
py -3.10 run.py figure1
py -3.10 plot.py 1
```

`run.py` builds its executable in the requested work directory. It does not need a pre-existing `build` folder. Model names are `nonheated`, `gamma500`, `gamma50`, `gamma10`, `delta1`, `delta2`, `delta3`. For each atmosphere, first iterate at 3 angles and then calculate the final field at 11 angles using the newly computed state:

```powershell
py -3.10 run.py nonheated --angles 3
py -3.10 run.py nonheated --angles 11 --initial work/nonheated_a3
py -3.10 run.py delta1 --angles 3
py -3.10 run.py delta1 --angles 11 --initial work/delta1_a3
py -3.10 plot.py 8 --delta 1
```

Repeat for `delta2` and `delta3`, then `plot.py 8` generates `1.pdf`, `2.pdf`, `3.pdf`. Compute the nonheated, gamma500, gamma50 and gamma10 models at 11 angles before running `plot.py 3`. Default new files go in `code/work`; use `--work` to put them elsewhere. Inspect commands without calculating with `run.py MODEL --dry-run`. Use `--output` for a new directory when continuing a previous state with `--initial`; existing directories are never overwritten. Atmosphere iterations may take hours. Exit code 3 means equilibrium was not reached. Plotting requires `status.json` to report convergence.

## Calculation scope

Deleting results also deletes the previously computed warm starts. The runner now starts from a grey temperature guess and a Planck field unless `--initial` is supplied. This is a fresh initialization, not an exact replay of the previous checkpoint sequence. Fresh heated convergence has not been established by the packaging check; convergence and agreement with the paper must be verified after solving. The configured iterations are limits, not a convergence guarantee. The retained code does not establish reproduction of every figure in the paper.

The current atmosphere implementation uses a unit free-free Gaunt factor. Recorded calculations used midpoint Compton energy integration (the default); `--compton moment_v1` selects the implemented energy-bin correction. `gamma10` retains its existing HAT heating projection and finite-difference temperature response options. These are implementation choices, not independently established author settings.

The optional `--lapack-dll PATH` selects the Windows LP64 OpenBLAS least-squares backend used in the latest power-law runs. No DLL is bundled; without this option the solver uses its built-in Jacobi backend. This can alter runtime and numerical roundoff. Do not use an ILP64 DLL.

Figure 8 plots interpolate in log energy, normalize by a cubic angular extension to mu=1, and draw curves only from the smallest computed positive angle node to 1. They do not extrapolate to mu=0. No paper curves are used as solver inputs or bundled overlays.
