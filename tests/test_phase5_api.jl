using Test
using LinearAlgebra

if isdefined(Main, :JuED)
    const EDMod = Main.JuED.EDMod
else
    include("../src/EDMain.jl")
    const EDMod = Main.EDMod
end

@testset "Explicit model constructors choose concrete types" begin
    onebody4 = zeros(Float64, 4, 4)
    two4 = zeros(Float64, 4, 4, 4, 4)
    two6 = zeros(Float64, 2, 2, 2, 2, 2, 2)

    @test EDMod.SpinlessListModel(2, 2, 2, onebody4, two4) isa EDMod.ModelTypesMod.ModelParams2DSpinlessList
    @test EDMod.SpinfulListModel(1, 1, 2, 2, onebody4, two4) isa EDMod.ModelTypesMod.ModelParams2DSpinList
    @test EDMod.SpinlessMomentumModel(2, 2, 2, onebody4, two6) isa EDMod.ModelTypesMod.ModelParams2DSpinless
    @test EDMod.SpinfulMomentumModel(1, 1, 2, 2, onebody4, two6) isa EDMod.ModelTypesMod.ModelParams2DSpin
    @test EDMod.TwoBandModel(2, 2, 1, zeros(Float64, 4, 4), two4) isa EDMod.ModelTypesMod.ModelParams2DTwoBand
end

@testset "Reusable sector solver matches legacy spinless list wrappers" begin
    onebody = Diagonal([1.0, 2.0, 3.0, 4.0]) |> Matrix
    twobody = zeros(Float64, 4, 4, 4, 4)
    model = EDMod.SpinlessListModel(2, 2, 2, onebody, twobody)

    hilbertspace = EDMod.BuildSector(model, 1)
    operator, dim = EDMod.BuildOperator(model, hilbertspace)
    config = EDMod.SolverConfig(1; return_vectors=1)
    result = EDMod.SolveSector(model, 1, config)
    legacy_vals, legacy_vecs = EDMod.DiagonalizeOneMomentum(model, 1, 1)

    @test dim == result.dim == length(hilbertspace.hilbert)
    @test size(operator, 1) == dim
    @test isapprox(result.values[1], legacy_vals[1]; atol=1e-8)
    @test isapprox(abs(dot(result.vectors[1], legacy_vecs[1])), 1.0; atol=1e-8)
end

@testset "SolveAllSectors matches DiagonalizeAllMomentum and supports matrixfree" begin
    onebody = zeros(Float64, 4, 4)
    twobody = zeros(Float64, 4, 4, 4, 4)
    twobody[1, 2, 2, 1] = 0.25
    model = EDMod.SpinlessListModel(2, 2, 2, onebody, twobody)

    results = EDMod.SolveAllSectors(model, EDMod.SolverConfig(1; return_vectors=0))
    legacy = EDMod.DiagonalizeAllMomentum(model, 1)
    matrixfree_result = EDMod.SolveSector(model, 1, EDMod.SolverConfig(1; return_vectors=1, matrixfree=true))
    sparse_result = EDMod.SolveSector(model, 1, EDMod.SolverConfig(1; return_vectors=1, matrixfree=false))
    expected = fill(NaN, 1, length(results))
    for (col, result) in enumerate(results)
        expected[1:length(result.values), col] = result.values
    end

    @test isequal(expected, legacy)
    @test isapprox(matrixfree_result.values[1], sparse_result.values[1]; atol=1e-8)
end

@testset "Two-band solver pipeline supports save-style vector collection" begin
    onebody = zeros(Float64, 4, 4)
    twobody = zeros(Float64, 4, 4, 4, 4)
    model = EDMod.TwoBandModel(2, 2, 1, onebody, twobody)

    elist, mocoeffs = EDMod.DiagonalizeAllMomentum(model, 1; numer_of_vectors=1, save=true)
    results = EDMod.SolveAllSectors(model, EDMod.SolverConfig(1; return_vectors=1))

    @test size(elist, 2) == 2
    @test Set(keys(mocoeffs)) == Set(result.momentum for result in results)
end

@testset "Unsupported spinful 6-index pipeline errors clearly" begin
    onebody = zeros(Float64, 4, 4)
    two6 = zeros(Float64, 2, 2, 2, 2, 2, 2)
    model = EDMod.SpinfulMomentumModel(1, 1, 2, 2, onebody, two6)
    hilbertspace = EDMod.BuildSector(model, 0)
    @test_throws ArgumentError EDMod.BuildOperator(model, hilbertspace)
end
