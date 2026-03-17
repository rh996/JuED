# JuED Refactor Program

## Purpose

This document turns the current codebase shortcomings into an implementation plan.
The goal is to make JuED easier to maintain, safer to extend, and faster on the
large exact-diagonalization workloads it targets.

## Main Problems Observed

1. Module structure is inconsistent.
   The entry module includes both `src/MomentumHilbertSpace2D_cache.jl` and
   `src/MomentumHilbertSpace2D.jl`, and both define `MomentumHilbertSpace2DMod`.
   This creates shadowing risk and makes it unclear which implementation is
   active.

2. The codebase is not organized as a standard Julia package.
   There is no visible `Project.toml`, and the tests are individual scripts
   rather than a conventional package test suite with shared fixtures and a
   single entry point.

3. The same algorithms are reimplemented many times.
   Momentum helpers, DFS basis builders, spin/two-band combiners, and
   diagonalization loops are duplicated across multiple files with minor
   variations.

4. Type usage is inconsistent and sometimes unsafe.
   The code mixes `Int32`, `Int64`, `UInt32`, `UInt64`, and hard-coded `Int64`
   state containers. Some arrays are built from `[]`, which can introduce `Any`
   or weak inference. In `Hamiltonian_list_constructor`, `row` is allocated with
   a dynamic pointer type but later filled via `Int32(...)`.

5. Hamiltonian assembly is correct in spirit but inefficient.
   Sparse construction repeatedly allocates temporary fermion objects, performs
   many dictionary lookups, and duplicates nearly the same two-pass CSC logic
   for different model types.

6. Progress reporting is threaded into hot loops.
   `ProgressMeter.next!` is called from inside threaded loops, which is not part
   of the scientific algorithm and can distort performance or complicate
   parallel execution.

7. RDM routines are the largest maintainability and scalability risk.
   `RDM1`, `RDM2`, and especially `RDM3` rebuild Hilbert spaces internally,
   allocate very large dense tensors, enumerate large operator index sets in
   ad hoc ways, and repeat the same operator-application pattern with slight
   variations.

8. Algorithm selection is coupled to data layout.
   The distinction between 4-index sparse interaction lists and 6-index momentum
   tensors is encoded by overloaded entry points rather than by explicit
   operator/model abstractions.

9. The public API is thin and the internal boundaries are weak.
   Entry-point functions directly choose index width, build Hilbert spaces,
   assemble sparse matrices, and invoke eigensolvers. That makes the code harder
   to profile, test, and reuse.

10. Performance validation infrastructure is missing.
    The repository has useful exploratory tests, but no benchmark suite, no
    profiling workflow, and no acceptance thresholds for large-sector runs.

## Refactor Goals

1. Establish one clear package/module layout.
2. Remove duplicate implementations and centralize shared algorithms.
3. Make basis generation, operator application, Hamiltonian assembly, and RDM
   evaluation explicit subsystems with stable interfaces.
4. Improve type stability and remove avoidable allocations.
5. Preserve current physics results while making large runs cheaper.
6. Add tests and benchmarks that detect correctness regressions and performance
   regressions.

## Program Overview

The work should be done in phases. Earlier phases reduce structural risk.
Later phases optimize the hot paths after the architecture is stable.

## Progress Update

Completed in the current branch:

1. Resolved the duplicate 2D momentum module path.
   `MomentumHilbertSpace2D_cache.jl` was removed and its memoized DFS logic was
   folded into the canonical `MomentumHilbertSpace2D.jl`.
2. Introduced a shared momentum utility module.
   `MomentumUtils.jl` now provides the canonical 1D and 2D momentum addition and
   subtraction helpers, and the basis modules import/re-export those functions.
3. Introduced a shared basis-builder module.
   `BasisBuilders.jl` now owns the generic particle-number DFS and the generic
   momentum-constrained DFS used by the 1D, 2D, spinful, and two-band Hilbert
   spaces.
4. Tightened several weakly typed Hilbert-space accumulators.
   The spinful and two-band basis combiners now initialize typed vectors instead
   of starting from bare `[]`.
5. Introduced explicit state and sparse-index helper utilities.
   `IndexTypes.jl` now centralizes state-width selection and sparse pointer-width
   selection, replacing repeated `if ... > 31` logic in the entry points and
   removing the unsigned pointer split in sparse Hamiltonian assembly.
6. Exposed explicit cache policy in the shared basis API.
   The shared basis builders and the public Hilbert-space builders now support
   `use_cache=true/false`, so memoization is no longer an implicit hidden choice.
7. Converted the remaining Phase 2 RDM helper collections to typed structures.
   Pair lists, triple lists, momentum-conserving index sets, and the `RDM2`
   cache dictionaries now use explicit tuple and pointer types instead of ad hoc
   dynamic containers.
8. Added focused structured regression tests for exact basis contents and cache
   policy behavior.
   The refactor test suite now checks cached vs uncached basis equality and
   exact small-basis contents for representative spinless and spinful cases.
9. Added a constructor-level verification pass after the type cleanup.
   `tests/test_list_constructor.jl` was rerun after the sparse pointer changes
   to confirm the Hamiltonian list-construction path still works.
10. Added a focused structured regression test.
   `tests/test_refactor_hilbert_spaces.jl` verifies the shared momentum helpers
   and checks that the refactored sector builders still partition the spinless,
   spinful, and two-band Hilbert spaces correctly on small systems.
11. Repaired the standalone 1D Hilbert-space test harness.
   `tests/test_hilbert_k1d.jl` now loads through `EDMain.jl`, so its dependency
   chain matches the actual module layout.

Still pending from the early phases:

1. Standard package metadata and a single canonical test runner.
2. Public/internal API separation and Hamiltonian-construction unification.
3. RDM architecture cleanup beyond collection typing, especially eliminating
   repeated Hilbert-space rebuilds and consolidating operator-application logic.

## Phase 0: Package and Repository Hygiene

1. Add a standard `Project.toml` and package metadata.
2. Create a single test entry point under Julia's standard `test/` layout, or
   standardize the existing `tests/` folder with a master runner.
3. Define development dependencies explicitly, including `KrylovKit`,
   `ProgressMeter`, `JLD2`, and any benchmarking tools.
4. Add a short contributor note that explains the expected package workflow,
   test commands, and performance-check workflow.

### Deliverables

- Standard package metadata
- Reproducible environment setup
- One command to run the whole test suite

## Phase 1: Module Boundary Cleanup

1. Resolve the duplicate `MomentumHilbertSpace2DMod` definitions. Completed.
   Choose one implementation path:
   - keep a cached implementation as the default and delete the duplicate file,
     or
   - rename the cached variant to a distinct module and make the choice explicit
     in the API.
2. Consolidate repeated momentum utility functions into one shared utility
   module. Completed for the 1D/2D momentum helpers.
3. Move all basis-space types into a single namespace with consistent exports.
4. Separate public entry points from internal implementation modules.

### Deliverables

- No duplicate module names
- One source of truth for momentum arithmetic
- Clear internal vs public API boundaries

## Phase 2: Basis Generation Unification

1. Introduce a shared basis-generation engine parameterized by: Completed.
   - orbital count
   - particle count
   - momentum constraint
   - bit layout strategy
2. Replace the repeated DFS implementations in: Completed for the current
   Hilbert-space files.
   - `HilbertSpace.jl`
   - `MomentumHilbertSpace1D.jl`
   - `MomentumHilbertSpace2D.jl`
   - `SpinHilbertSpace.jl`
   - `SpinMomentumHilbertSpace1D.jl`
   - `SpinMomentumHilbertSpace2D.jl`
   - `TwoBandMomentumHilbertSpace2D.jl`
3. Make caching a configurable feature of the basis builder instead of a
   separate shadow implementation. Completed.
   The shadow implementation was removed, memoization lives in the shared basis
   builder, and `use_cache=true/false` is now exposed in the shared/public basis
   builder APIs.
4. Introduce a dedicated state-index type alias and enforce the same type across
   Hilbert spaces, dictionaries, and sparse rows. Completed for the current
   basis/setup layer.
   State-width and pointer-width selection go through shared helpers, `ToDict`
   can build typed index maps, the sparse constructors use shared pointer
   helpers, and the Phase 2 RDM/cache collections now use explicit tuple/index
   types.

### Deliverables

- One generalized DFS basis builder
- Explicit cache policy
- Consistent state/index typing

### Phase 2 Follow-Up

Phase 2 is complete for the basis-generation and sector-typing scope.
The remaining follow-up items now belong to later phases:

1. Further RDM redesign belongs to Phase 6.
2. Broader sparse-assembly cleanup belongs to Phase 4.
3. Package/test-runner standardization belongs to Phase 0.

## Phase 3: Fermion Operator Kernel Cleanup

1. Generalize `FermionOperator` over the state integer type instead of hardcoding
   `Int64`.
2. Separate pure operator kernels from mutable convenience wrappers.
3. Add micro-tests for creation/annihilation parity, occupancy checks, and
   corner cases near the maximum supported bit width.
4. Document the basis ordering convention and reverse-index mapping currently
   embedded throughout the code.

### Deliverables

- Typed fermion operator kernel
- Explicit basis-ordering documentation
- Unit tests for sign conventions

## Phase 4: Hamiltonian Construction Refactor

1. Extract a generic sparse-column assembly engine with a model-specific callback
   for generating connected states.
2. Remove duplicated two-pass CSC construction logic between the 4-index and
   6-index interaction paths.
3. Precompute reusable operator index maps for momentum-conserving scattering
   channels.
4. Remove `ProgressMeter` calls from hot threaded loops or gate them behind a
   debug/performance flag.
5. Audit integer conversions so `row`, `indptr`, and dictionary indices use one
   consistent type.
6. Add a matrix-free Hamiltonian application path as a first-class API alongside
   explicit sparse construction, based on the pattern already explored in the
   vector-product tests.

### Deliverables

- One Hamiltonian assembly framework
- Optional matrix-free action for large sectors
- Reduced allocation pressure in hot loops

## Phase 5: Public API Simplification

1. Replace overload-heavy `InputModel` construction with explicit model types or
   constructors that make the interaction representation obvious.
2. Split "build basis", "build Hamiltonian", and "solve" into separately callable
   steps.
3. Introduce a solver configuration object for tolerances, number of eigenpairs,
   return-vectors behavior, and sparse vs matrix-free mode.
4. Make sector iteration reusable for all model families rather than duplicating
   `DiagonalizeAllMomentum` loops.

### Deliverables

- Smaller, more explicit API surface
- Reusable sector solver pipeline
- Clear configuration semantics

## Phase 6: Density Matrix Redesign

1. Stop rebuilding Hilbert spaces inside every RDM call.
   Accept a prepared sector object or cached workspace where possible.
2. Factor out the shared operator-application pattern used by `RDM1`, `RDM2`,
   `RDM3`, and cache generation.
3. Replace ad hoc tuple-generation code with reusable momentum-filter utilities.
4. Make caching a structured subsystem with explicit filenames, schema versioning,
   and invalidation rules.
5. Re-evaluate dense tensor outputs for `RDM2` and `RDM3`.
   For larger systems, expose sparse or symmetry-reduced representations to avoid
   allocating full `norb^4` and `norb^6` tensors when they are not needed.
6. Add correctness tests that compare optimized RDM code to naive reference
   implementations on small systems.

### Deliverables

- Shared RDM engine
- Reusable cache format
- Lower memory footprint for higher-body density matrices

## Phase 7: Performance and Memory Work

1. Add `BenchmarkTools` benchmarks for:
   - Hilbert-space generation
   - Hamiltonian assembly
   - sparse eigensolve
   - matrix-free eigensolve
   - `RDM1`, `RDM2`, and `RDM3`
2. Profile allocations in the operator kernels and sparse assembly loops.
3. Preallocate reusable thread-local workspaces for repeated operator
   applications.
4. Replace repeated `Dict` lookups with denser indexing structures where sector
   ordering makes that possible.
5. Define target problem sizes and record expected runtime and memory envelopes.

### Deliverables

- Benchmark suite
- Allocation profile baselines
- Performance targets for future changes

## Phase 8: Testing and Verification

1. Convert the existing exploratory scripts into structured tests with assertions.
2. Add cross-checks between:
   - sparse Hamiltonian vs matrix-free action
   - cached vs uncached basis generation
   - optimized RDMs vs naive RDMs
   - old API outputs vs refactored API outputs on small systems
3. Add regression tests for hermiticity, particle-number conservation, and
   momentum-sector consistency.
4. Add CI for the test suite and, if runtime permits, a small benchmark smoke
   test.

### Deliverables

- Reproducible correctness coverage
- Regression protection for future optimization work

## Recommended Execution Order

1. Phase 0
2. Phase 1
3. Phase 2
4. Phase 3
5. Phase 4
6. Phase 5
7. Phase 6
8. Phase 8
9. Phase 7

Performance work should start only after the architectural duplication is
removed, otherwise profiling results will be noisy and the same optimization
will have to be repeated in several places.

## Immediate First Refactors

These are the highest-leverage changes to start with:

1. Remove the duplicate 2D momentum module definitions. Completed.
2. Standardize the package structure and test entry point. In progress.
3. Unify the DFS basis builders behind one generic implementation.
4. Normalize integer and container types across basis states and sparse data.
5. Extract shared Hamiltonian assembly scaffolding.
6. Redesign RDM code around a cached sector workspace.

## Success Criteria

The refactor program is complete when:

1. The project has one clear package structure and one canonical module graph.
2. Shared basis-generation and Hamiltonian code paths replace the existing
   duplicated variants.
3. Public APIs express models, sectors, and solver settings explicitly.
4. Small-system results are unchanged relative to the current implementation.
5. Large-system runs allocate less memory and complete faster in the main hot
   paths.
6. Benchmarks and tests catch both correctness and performance regressions.
