# JuED

JuED is a Julia exact-diagonalization package for fermionic lattice models with
momentum-, spin-, and two-band sector support.

The maintained package surface covers:

- symmetry-reduced Hilbert-space construction
- sparse and matrix-free Hamiltonian application
- momentum-sector eigensolving with `KrylovKit`
- one-, two-, and three-body reduced density matrices

## Installation

Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/rh996/JuED.git")
```

Then load the package:

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

## Main Concepts

### Basis spaces

Use `BasisSpaces` for the maintained Hilbert-space namespace:

- `BasisSpaces.GeneralHilbertSpace`
- `BasisSpaces.MomentumHilbertSpace1D`
- `BasisSpaces.MomentumHilbertSpace2D`
- `BasisSpaces.SpinHilbertSpace`
- `BasisSpaces.SpinMomentumHilbertSpace1D`
- `BasisSpaces.SpinMomentumHilbertSpace2D`
- `BasisSpaces.TwoBandMomentumHilbertSpace2D`

Each basis object is built with `BasisSpaces.build_hilbert!(...)`.

### Models

The preferred model constructors are:

- `SpinlessListModel`
- `SpinlessMomentumModel`
- `SpinfulListModel`
- `SpinfulMomentumModel`
- `TwoBandModel`

### Solver pipeline

The reusable solve pipeline is:

1. `BuildSector(model, momentum)`
2. `BuildOperator(model, hilbertspace)`
3. `SolveSector(model, momentum, config)`
4. `SolveAllSectors(model, config)`

### Reduced density matrices

For repeated work in the same momentum sector, keep an `RDMWorkspace`:

```julia
workspace = RDMWorkspace(model, 0)
rdm3_compact = RDM3Compact(workspace, ground_vector)
```

Compact higher-body RDMs can be persisted:

```julia
workspace = RDMWorkspace(model, 0)
file = compact_rdm_filename(workspace, 3)
RDM3Compact(workspace, ground_vector; savefile=file)
rdm3_loaded = load_compact_rdm3(workspace; file=file)
rdm3_dense = todense(rdm3_loaded)
```

## Package Layout

- `src/JuED.jl`: package entrypoint
- `src/PublicAPI.jl`: package-facing API
- `src/Internal.jl`: subsystem include graph
- `src/EDMain.jl`: compatibility facade for older `EDMod`-based code paths
- `test/runtests.jl`: canonical package test suite

## Local Docs Build

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```
