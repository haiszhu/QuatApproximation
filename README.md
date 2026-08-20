# QuatApproximation (outdated)

Quaternion-approximation reference code (Fortran + MATLAB MEX on Mac).

## Context
This code was written by Codex and Claude Code as a testbed to reproduce parts of:
- Hai Zhu and Shravan Veerapaneni, *High-order close evaluation of Laplace layer potentials: a differential geometric approach*, SIAM Journal on Scientific Computing, 44(3): A1381–A1405 (2022).
  https://doi.org/10.1137/21M1408331
- Shidong Jiang and Hai Zhu, *Recursive reduction quadrature for the evaluation of Laplace layer potentials in three dimensions*, arXiv:2411.08342 (2024).
  https://arxiv.org/abs/2411.08342

One purpose is to provide the quaternion-approximation and solid-spherical-harmonic building blocks (`q^{(n,m)}_i`, `omega^{(n,m)}_i`, harmonic basis values + gradients, Koornwinder/Vioreanu and tensor-product reference-element geometry, ...) in r64 and r128 so they can be combined with [LineQuaaadrature](https://github.com/haiszhu/LineQuaaadrature) to serve the close-evaluation chain in [qotential](https://github.com/haiszhu/qotential).

## Layout
- `src/` Fortran source modules
- `matlab/` mwrap interface (`QuatApproximation.mw`) and generated MEX wrappers
- `test/` Fortran test programs
- `build/` compiled artifacts

## Build
From this folder:

```bash
make
```

This builds:
- static library: `build/libQuatApproximation.a`
- MATLAB MEX: `matlab/QuatApproximation_mex.mexmaca64` (on Apple Silicon macOS)

### Windows

MATLAB R2024a or later is required for MinGW Fortran MEX support. In
MATLAB's Add-On Explorer, install **MATLAB Support for MinGW-w64
C/C++/Fortran Compiler**, then configure both compilers:

```matlab
mex -setup C
mex -setup FORTRAN
```

The HDF5-enabled build also requires an HDF5 development installation
containing `include/hdf5.h`, `lib/libhdf5.dll.a`, and
`bin/libhdf5.dll`. Build and test from PowerShell with:

```powershell
make -f makefile.windows
make -f makefile.windows test
```

The Windows build stages the HDF5 runtime DLLs beside the MEX file so it
can be loaded from an ordinary MATLAB session without changing `PATH`.

Defaults target MATLAB R2024b, its bundled MinGW toolchain, and the HDF5
development files bundled with the local Octave installation. Override
any location when needed:

```powershell
make -f makefile.windows MATLABROOT="D:/MATLAB/R2024b" `
    MINGWROOT="D:/MathWorks/MinGW" HDF5_ROOT="D:/Libraries/hdf5"
```

Remove generated Windows products with:

```powershell
make -f makefile.windows clean
```

## Run Fortran Tests

### Paraboloid approximation conditioning
Probes the condition number of the harmonic-approximation matrix on a paraboloid patch across mean-curvature values:
```bash
make test_paraboloid_cond
./build/test_paraboloid_cond
```

### Paraboloid refinement (Fig 5.1 right reproducer)
Builds the 4×4 quaternion-approximation block matrix, solves at increasing panel orders + refinement levels, and emits the convergence + condition-number tables:
```bash
make test_paraboloid_refinement
./build/test_paraboloid_refinement
```
PDFs produced via gnuplot: `paraboloid_approx_error_vs_refinement.pdf`, `paraboloid_cond_num_vs_refinement.pdf`.

Optional thread control:
```bash
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 ./build/test_paraboloid_refinement
```

## Clean
```bash
make clean
```

## Roadmap
Open items are tracked in [`PLAN.md`](PLAN.md).  Notably, the local-to-local translation operator (transferring a solid spherical harmonic expansion from one center to another) is not yet implemented; it is required for refinement-based close-evaluation strategies and is the next major item.

## Notes
- Compiler settings are in `Makefile` (currently `gfortran-15`, `gcc-15`).
- MEX wrappers are generated from `matlab/QuatApproximation.mw` using `mwrap`.
- `src/SUBROUTINES.md` is the authoritative per-routine inventory (signatures, purposes, r64/r128 status).

## To do list

* Slm basis, locloc, etc (redesign with both dense mat vec and sparse for loop)
