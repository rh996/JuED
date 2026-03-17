# Benchmark Targets

This directory defines the maintained JuED benchmark workloads and the expected
performance envelope for each one. These are soft targets for developer use,
not CI gates.

## Canonical Command

```bash
julia --project=. benchmarks/runbenchmarks.jl
```

The benchmark runner writes:

- `benchmarks/results/latest.md`
- `benchmarks/results/latest.tsv`

## Workloads

- `kernels/apply_operator_string`
  Small fermion-kernel microbenchmark used to watch allocation regressions in
  the typed operator path.
- `basis/general_hilbert_8o_3p`
  Generic particle-number Hilbert-space generation.
- `basis/spinless_k2d_2x2_n2_k1`
  Symmetry-reduced 2D spinless momentum-sector basis generation.
- `basis/spinful_k1d_4o_n2n2_k0`
  1D spinful momentum-sector basis generation.
- `assembly/list_spinless_2x2_n2`
  Sparse CSC assembly for the 4-index interaction path.
- `assembly/momentum_spinless_2x2_n2`
  Sparse CSC assembly for the 6-index momentum interaction path.
- `solve/sparse_spinless_list_2x2_n2`
  Sparse eigensolve on a prebuilt operator.
- `solve/matrixfree_spinless_list_2x2_n2`
  Matrix-free eigensolve on the equivalent sector/problem.
- `rdm/rdm1_spinless_2x2_n2`
  `RDM1` workspace evaluation.
- `rdm/rdm2compact_spinless_2x2_n2`
  Compact `RDM2` workspace evaluation.
- `rdm/rdm3compact_spinless_5x1_n3`
  Compact `RDM3` workspace evaluation on a small nontrivial sector.

## Envelope Policy

When you compare new results to the latest baseline:

- Kernel and basis workloads should remain allocation-light.
- Assembly workloads should not regress by more than one order of magnitude in
  median time or allocated bytes without a clear algorithmic reason.
- Sparse and matrix-free solves should be compared to each other on the same
  workload rather than treated as interchangeable absolute targets.
- `RDM2` and `RDM3` compact paths should remain the primary benchmarked higher-
  body workloads; dense reconstruction is intentionally excluded from the
  maintained benchmark suite because it scales with tensor materialization.

## Current Soft Envelopes

The following envelopes were recorded from the maintained benchmark runner on
the current branch with one Julia thread and the default benchmark settings.

| Benchmark | Observed Median (ns) | Observed Memory (bytes) | Soft Ceiling (ns) | Soft Ceiling (bytes) |
| --- | ---: | ---: | ---: | ---: |
| `kernels/apply_operator_string` | 0 | 0 | 100 | 0 |
| `basis/general_hilbert_8o_3p` | 1396 | 5776 | 10000 | 8192 |
| `basis/spinless_k2d_2x2_n2_k1` | 500 | 1552 | 5000 | 4096 |
| `basis/spinful_k1d_4o_n2n2_k0` | 3021 | 12688 | 15000 | 16384 |
| `assembly/list_spinless_2x2_n2` | 31084 | 4256 | 100000 | 8192 |
| `assembly/momentum_spinless_2x2_n2` | 28750 | 3520 | 100000 | 8192 |
| `solve/sparse_spinless_list_2x2_n2` | 3896 | 21248 | 20000 | 32768 |
| `solve/matrixfree_spinless_list_2x2_n2` | 40958 | 25744 | 100000 | 32768 |
| `rdm/rdm1_spinless_2x2_n2` | 125 | 336 | 2000 | 1024 |
| `rdm/rdm2compact_spinless_2x2_n2` | 14146 | 1168 | 50000 | 4096 |
| `rdm/rdm3compact_spinless_5x1_n3` | 12667 | 1360 | 50000 | 4096 |
