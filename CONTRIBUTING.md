# Contributing

## Environment

Use the repository project directly:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Canonical Test Command

Run the maintained package regression suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This executes `test/runtests.jl`. The scripts under `tests/` are still useful
for local experiments, but they are not all part of the maintained package
suite.

## Benchmark Command

Run the maintained benchmark harness with:

```bash
julia --project=. benchmarks/runbenchmarks.jl
```

This writes a Markdown summary and TSV table under `benchmarks/results/`. Use
`benchmarks/TARGETS.md` as the reference workload list and soft-envelope guide
when comparing runs.

## Development Notes

- Prefer updating the maintained regression tests when changing public behavior.
- Keep new package-facing APIs in `src/EDMain.jl` or `src/JuED.jl`.
- Keep lower-level implementation changes in the subsystem files under `src/`.
- When adding a new dependency, update `Project.toml`.
- When adding a new maintained regression, wire it into `test/runtests.jl`.
