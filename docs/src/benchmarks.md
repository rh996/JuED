# Benchmarks

JuED includes a maintained benchmark harness based on `BenchmarkTools`.

## Run Benchmarks

```bash
julia --project=. benchmarks/runbenchmarks.jl
```

The runner writes:

- `benchmarks/results/latest.md`
- `benchmarks/results/latest.tsv`

## Covered Workloads

The maintained suite benchmarks:

- typed fermion-kernel application
- representative basis generation workloads
- sparse Hamiltonian assembly
- sparse and matrix-free eigensolve paths
- compact `RDM1`, `RDM2`, and `RDM3` workloads

## Targets

The benchmark workload list and current soft runtime/memory envelopes live in:

- `benchmarks/TARGETS.md`

Those targets are not CI gates. They are intended to make optimization work and
regression checking reproducible across local runs.
