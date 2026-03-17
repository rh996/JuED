# JuED

JuED is a Julia exact-diagonalization package for fermionic lattice models. It
provides:

- symmetry-reduced Hilbert-space construction
- sparse and matrix-free Hamiltonian application
- momentum-sector eigensolvers built on `KrylovKit`
- one-, two-, and three-body reduced density matrices

## Installation

At the moment, install JuED directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/rh996/JuED.git")
```

Then load it with:

```julia
using JuED
```

## Quick Start

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
rdm2 = RDM2(model, ground_vector, 0; representation=:compact)
```

For repeated work in the same momentum sector, reuse a workspace:

```julia
workspace = RDMWorkspace(model, 0)
rdm3_compact = RDM3Compact(workspace, ground_vector)
```

If you want to persist compact higher-body RDMs:

```julia
workspace = RDMWorkspace(model, 0)
file = compact_rdm_filename(workspace, 3)
RDM3Compact(workspace, ground_vector; savefile=file)
rdm3_loaded = load_compact_rdm3(workspace; file=file)
rdm3_dense = todense(rdm3_loaded)
```

## Public API

The preferred package surface is:

- `BasisSpaces` for Hilbert-space types and basis builders
- `SpinlessListModel`, `SpinlessMomentumModel`, `SpinfulListModel`,
  `SpinfulMomentumModel`, `TwoBandModel` for explicit model construction
- `BuildSector`, `BuildOperator`, `SolveSector`, `SolveAllSectors` for the
  reusable solve pipeline
- `SolverConfig` for eigensolver settings
- `RDMWorkspace`, `RDM1`, `RDM2`, `RDM3`, `RDM2Compact`, `RDM3Compact` for
  reduced-density-matrix workflows

Backward-compatible wrappers such as `InputModel`,
`DiagonalizeOneMomentum`, and `DiagonalizeAllMomentum` are still available.

## Development

For local development:

```bash
git clone git@github.com:rh996/JuED.git
cd JuED
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the maintained test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run the maintained benchmark harness with:

```bash
julia --project=. benchmarks/runbenchmarks.jl
```

See `CONTRIBUTING.md` for the development workflow and
`benchmarks/TARGETS.md` for the maintained benchmark workloads.

Build the local documentation site with:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## License

JuED is distributed under the MIT License. See `LICENSE`.
