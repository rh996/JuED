# JuED

JuED is a Julia exact-diagonalization codebase for fermionic lattice models with
momentum-, spin-, and two-band sector support. The maintained package surface
currently covers:

- symmetry-reduced Hilbert-space construction
- sparse and matrix-free Hamiltonian application
- momentum-sector eigensolving with `KrylovKit`
- one-, two-, and three-body reduced density matrices

## Install

```bash
git clone <repo-url>
cd JuED
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The repository now includes a standard `Project.toml`, package entrypoint at
`src/JuED.jl`, and canonical package test runner at `test/runtests.jl`.

## Basic Usage

```julia
using JuED

nparticle = 2
Nkx = 2
Nky = 2
onebody = zeros(Float64, 4, 4)
twobody = zeros(Float64, 4, 4, 4, 4)

model = SpinlessListModel(nparticle, Nkx, Nky, onebody, twobody)
config = SolverConfig(1; return_vectors=1)

result = SolveSector(model, 0, config)
ground_energy = result.values[1]
ground_vector = result.vectors[1]

rdm1 = RDM1(model, ground_vector, 0)
rdm2_compact = RDM2(model, ground_vector, 0; representation=:compact)
```

Backward-compatible wrappers such as `InputModel`, `DiagonalizeOneMomentum`, and
`DiagonalizeAllMomentum` are still available, but the preferred public API is:

- explicit model constructors:
  `SpinlessListModel`, `SpinlessMomentumModel`, `SpinfulListModel`,
  `SpinfulMomentumModel`, `TwoBandModel`
- reusable solver pipeline:
  `BuildSector`, `BuildOperator`, `SolveSector`, `SolveAllSectors`
- configurable solver settings:
  `SolverConfig`

## Test Workflow

Run the maintained package regression suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This executes `test/runtests.jl`, which covers the maintained regression tests
for:

- Hilbert-space refactors
- fermion-operator kernels
- Hamiltonian sparse vs matrix-free parity
- RDM workspace and compact paths
- public API regression

The `tests/` directory still contains older exploratory scripts and model-
specific experiments. Those are useful for local investigation, but they are
not the canonical package test suite.

## Repository Layout

- `src/JuED.jl`: package entrypoint
- `src/EDMain.jl`: public ED API and compatibility wrappers
- `src/HamiltonianConstructor.jl`: sparse and matrix-free Hamiltonian builders
- `src/DensityMatrices.jl`: RDM workspace and compact/dense RDM paths
- `test/runtests.jl`: canonical package test runner
- `tests/`: maintained regression files plus older exploratory scripts

## Development

See `CONTRIBUTING.md` for the expected development workflow.
