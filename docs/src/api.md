# API

This page summarizes the maintained package-facing API.

## Basis Namespace

Use `BasisSpaces` as the single maintained Hilbert-space namespace.

Common entry points:

- `BasisSpaces.build_hilbert!(hilbertspace; use_cache=true)`
- `BasisSpaces.state_index_map(hilbertspace, Int)`
- `BasisSpaces.momentum_add_1d`, `BasisSpaces.momentum_sub_1d`
- `BasisSpaces.momentum_add_2d`, `BasisSpaces.momentum_sub_2d`

## Model Constructors

Preferred constructors:

- `SpinlessListModel(nparticle, Nkx, Nky, onebody, twobody4)`
- `SpinlessMomentumModel(nparticle, Nkx, Nky, onebody, twobody6)`
- `SpinfulListModel(nalpha, nbeta, Nkx, Nky, onebody, twobody4)`
- `SpinfulMomentumModel(nalpha, nbeta, Nkx, Nky, onebody, twobody6)`
- `TwoBandModel(nparticle, Nkx, Nky, onebody, twobody4)`

Compatibility wrappers:

- `InputModel(...)`
- `InputTwoBandModel(...)`

## Solver API

The maintained solver entry points are:

- `SolverConfig(neigen; return_vectors=1, matrixfree=false, tol=1e-6, maxiter=1000, which=:SR, ishermitian=true)`
- `BuildSector(model, momentum; use_cache=true)`
- `BuildOperator(model, hilbertspace; matrixfree=false)`
- `SolveSector(model, momentum, config=SolverConfig(1); use_cache=true)`
- `SolveAllSectors(model, config=SolverConfig(1); use_cache=true)`

Legacy diagonalization wrappers still exist:

- `DiagonalizeOneMomentum(...)`
- `DiagonalizeAllMomentum(...)`

## RDM API

Workspace-oriented API:

- `RDMWorkspace(model, momentum; use_cache=true)`
- `RDM1(workspace, coeffs)`
- `RDM2(workspace, coeffs; representation=:dense | :compact)`
- `RDM3(workspace, coeffs; representation=:dense | :compact)`
- `RDM2Compact(workspace, coeffs)`
- `RDM3Compact(workspace, coeffs)`
- `todense(compact_rdm)`

Compact persistence helpers:

- `compact_rdm_filename(workspace, order; dir="./data")`
- `save_compact_rdm(compact_rdm, workspace; file=...)`
- `load_compact_rdm2(workspace; file=...)`
- `load_compact_rdm3(workspace; file=...)`

Transition-cache helper:

- `RDM2_cache(workspace; file=...)`

## Compatibility Layer

`JuED.EDMod` is still available for older code and tests, but it is a
compatibility facade. New code should prefer `using JuED` and the exported
package API directly.
