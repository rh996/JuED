using BenchmarkTools
using Dates
using JuED
using KrylovKit
using LinearAlgebra
using Printf
using Random

const BENCHMARK_SECONDS = parse(Float64, get(ENV, "JUED_BENCH_SECONDS", "0.75"))
const BENCHMARK_SAMPLES = parse(Int, get(ENV, "JUED_BENCH_SAMPLES", "10"))
const RESULT_DIR = joinpath(@__DIR__, "results")
const SUMMARY_DATE_FORMAT = dateformat"yyyy-mm-dd HH:MM:SS"

BenchmarkTools.DEFAULT_PARAMETERS.seconds = BENCHMARK_SECONDS
BenchmarkTools.DEFAULT_PARAMETERS.samples = BENCHMARK_SAMPLES

function normalized_state(dim::Int, seed::Int)
    rng = MersenneTwister(seed)
    coeffs = randn(rng, ComplexF64, dim)
    coeffs ./= norm(coeffs)
    return coeffs
end

function build_list_model()
    onebody = Diagonal([1.0, 2.0, 3.0, 4.0]) |> Matrix
    twobody = zeros(Float64, 4, 4, 4, 4)
    twobody[1, 2, 2, 1] = 0.25
    twobody[2, 1, 1, 2] = 0.25
    return JuED.SpinlessListModel(2, 2, 2, onebody, twobody), 1
end

function build_momentum_model()
    onebody = zeros(Float64, 4, 4)
    twobody = zeros(Float64, 2, 2, 2, 2, 2, 2)
    twobody[1, 1, 1, 1, 1, 1] = 0.5
    return JuED.SpinlessMomentumModel(2, 2, 2, onebody, twobody), 1
end

function build_rdm3_model()
    onebody = zeros(ComplexF64, 5, 5)
    twobody = zeros(ComplexF64, 5, 5, 5, 5)
    return JuED.SpinlessListModel(3, 5, 1, onebody, twobody)
end

function collect_rows(group::BenchmarkTools.BenchmarkGroup, prefix::Vector{String}=String[])
    rows = NamedTuple[]
    for key in sort!(collect(keys(group)); by=x -> string(x))
        value = group[key]
        if value isa BenchmarkTools.BenchmarkGroup
            append!(rows, collect_rows(value, [prefix; string(key)]))
        else
            estimate = median(value)
            push!(rows, (
                name=join([prefix; string(key)], "/"),
                time_ns=estimate.time,
                memory_bytes=estimate.memory,
                allocs=estimate.allocs,
            ))
        end
    end
    return rows
end

function write_summary(rows, output_dir::String)
    mkpath(output_dir)
    tsv_file = joinpath(output_dir, "latest.tsv")
    md_file = joinpath(output_dir, "latest.md")

    open(tsv_file, "w") do io
        println(io, "name\ttime_ns\tmemory_bytes\tallocs")
        for row in rows
            println(io, "$(row.name)\t$(row.time_ns)\t$(row.memory_bytes)\t$(row.allocs)")
        end
    end

    open(md_file, "w") do io
        println(io, "# JuED Benchmark Summary")
        println(io)
        println(io, "- Generated: $(Dates.format(now(), SUMMARY_DATE_FORMAT))")
        println(io, "- Julia threads: $(Threads.nthreads())")
        println(io, "- BenchmarkTools seconds/sample target: $(BENCHMARK_SECONDS)")
        println(io, "- BenchmarkTools samples: $(BENCHMARK_SAMPLES)")
        println(io)
        println(io, "| Benchmark | Median Time (ns) | Memory (bytes) | Allocs |")
        println(io, "| --- | ---: | ---: | ---: |")
        for row in rows
            println(io, @sprintf("| %s | %d | %d | %d |", row.name, row.time_ns, row.memory_bytes, row.allocs))
        end
    end

    return md_file, tsv_file
end

function build_suite()
    suite = BenchmarkGroup()

    list_model, list_momentum = build_list_model()
    momentum_model, momentum_sector = build_momentum_model()
    rdm3_model = build_rdm3_model()
    list_workspace = JuED.RDMWorkspace(list_model, list_momentum)
    list_coeffs = normalized_state(length(list_workspace.hilbert), 1234)
    rdm3_workspace = JuED.RDMWorkspace(rdm3_model, 0)
    rdm3_coeffs = normalized_state(length(rdm3_workspace.hilbert), 5678)

    suite["kernels"]["apply_operator_string"] = @benchmarkable JuED.EDMod.FermionOperatorMod.apply_operator_string(Int32(0x0003), (3, 4), (1, 2))

    suite["basis"]["general_hilbert_8o_3p"] = @benchmarkable begin
        hilbertspace = JuED.BasisSpaces.GeneralHilbertSpace{Int32}(3, 8, Int32[])
        JuED.BasisSpaces.build_hilbert!(hilbertspace; use_cache=true)
    end

    suite["basis"]["spinless_k2d_2x2_n2_k1"] = @benchmarkable begin
        hilbertspace = JuED.BasisSpaces.MomentumHilbertSpace2D{Int32}(2, 2, 2, 1, Int32[])
        JuED.BasisSpaces.build_hilbert!(hilbertspace; use_cache=true)
    end

    suite["basis"]["spinful_k1d_4o_n2n2_k0"] = @benchmarkable begin
        hilbertspace = JuED.BasisSpaces.SpinMomentumHilbertSpace1D{Int64}(2, 2, 4, 0, Int64[])
        JuED.BasisSpaces.build_hilbert!(hilbertspace; use_cache=true)
    end

    suite["assembly"]["list_spinless_2x2_n2"] = @benchmarkable begin
        JuED.BuildOperator($list_model, hilbertspace; matrixfree=false)
    end setup=(begin
        hilbertspace = JuED.BuildSector($list_model, $list_momentum)
    end)

    suite["assembly"]["momentum_spinless_2x2_n2"] = @benchmarkable begin
        JuED.BuildOperator($momentum_model, hilbertspace; matrixfree=false)
    end setup=(begin
        hilbertspace = JuED.BuildSector($momentum_model, $momentum_sector)
    end)

    suite["solve"]["sparse_spinless_list_2x2_n2"] = @benchmarkable begin
        KrylovKit.eigsolve(operator, dim, 1, :SR; maxiter=200, tol=1e-8, ishermitian=true)
    end setup=(begin
        hilbertspace = JuED.BuildSector($list_model, $list_momentum)
        operator, dim = JuED.BuildOperator($list_model, hilbertspace; matrixfree=false)
    end)

    suite["solve"]["matrixfree_spinless_list_2x2_n2"] = @benchmarkable begin
        KrylovKit.eigsolve(action, dim, 1, :SR; maxiter=200, tol=1e-8, ishermitian=true)
    end setup=(begin
        hilbertspace = JuED.BuildSector($list_model, $list_momentum)
        action, dim = JuED.BuildOperator($list_model, hilbertspace; matrixfree=true)
    end)

    suite["rdm"]["rdm1_spinless_2x2_n2"] = @benchmarkable JuED.RDM1($list_workspace, $list_coeffs)
    suite["rdm"]["rdm2compact_spinless_2x2_n2"] = @benchmarkable JuED.RDM2Compact($list_workspace, $list_coeffs)
    suite["rdm"]["rdm3compact_spinless_5x1_n3"] = @benchmarkable JuED.RDM3Compact($rdm3_workspace, $rdm3_coeffs)

    return suite
end

function main()
    suite = build_suite()
    results = run(suite; verbose=true)
    rows = collect_rows(results)
    md_file, tsv_file = write_summary(rows, RESULT_DIR)
    println("Wrote benchmark summary to $(md_file)")
    println("Wrote benchmark table to $(tsv_file)")
end

main()
