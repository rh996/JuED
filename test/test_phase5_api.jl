@testset "Compatibility constructors dispatch to the explicit model types" begin
    onebody4 = zeros(Float64, 4, 4)
    two4 = zeros(Float64, 4, 4, 4, 4)
    two6 = zeros(Float64, 2, 2, 2, 2, 2, 2)

    @test typeof(InputModel(2, 2, 2, onebody4, two4)) === typeof(SpinlessListModel(2, 2, 2, onebody4, two4))
    @test typeof(InputModel(1, 1, 2, 2, onebody4, two4)) === typeof(SpinfulListModel(1, 1, 2, 2, onebody4, two4))
    @test typeof(InputModel(2, 2, 2, onebody4, two6)) === typeof(SpinlessMomentumModel(2, 2, 2, onebody4, two6))
    @test typeof(InputModel(1, 1, 2, 2, onebody4, two6)) === typeof(SpinfulMomentumModel(1, 1, 2, 2, onebody4, two6))
    @test typeof(InputTwoBandModel(2, 2, 1, onebody4, two4)) === typeof(TwoBandModel(2, 2, 1, onebody4, two4))
end

@testset "Reusable sector solver matches legacy spinless list wrappers" begin
    onebody = Diagonal([1.0, 2.0, 3.0, 4.0]) |> Matrix
    twobody = zeros(Float64, 4, 4, 4, 4)
    model = SpinlessListModel(2, 2, 2, onebody, twobody)

    hilbertspace = BuildSector(model, 1)
    operator, dim = BuildOperator(model, hilbertspace)
    config = SolverConfig(1; return_vectors=1)
    result = SolveSector(model, 1, config)
    legacy_vals, legacy_vecs = DiagonalizeOneMomentum(model, 1, 1)

    @test dim == result.dim == length(hilbertspace.hilbert)
    @test size(operator, 1) == dim
    @test isapprox(result.values[1], legacy_vals[1]; atol=1e-8)
    @test isapprox(abs(dot(result.vectors[1], legacy_vecs[1])), 1.0; atol=1e-8)
end

@testset "SolveAllSectors matches DiagonalizeAllMomentum and supports matrixfree" begin
    onebody = zeros(Float64, 4, 4)
    twobody = zeros(Float64, 4, 4, 4, 4)
    twobody[1, 2, 2, 1] = 0.25
    model = SpinlessListModel(2, 2, 2, onebody, twobody)

    results = SolveAllSectors(model, SolverConfig(1; return_vectors=0))
    legacy = DiagonalizeAllMomentum(model, 1)
    matrixfree_result = SolveSector(model, 1, SolverConfig(1; return_vectors=1, matrixfree=true))
    sparse_result = SolveSector(model, 1, SolverConfig(1; return_vectors=1, matrixfree=false))
    expected = fill(NaN, 1, length(results))
    for (col, result) in enumerate(results)
        expected[1:length(result.values), col] = result.values
    end

    @test isapprox(expected, legacy; atol=1e-12, nans=true)
    @test isapprox(matrixfree_result.values[1], sparse_result.values[1]; atol=1e-8)
end

@testset "Two-band solver pipeline supports save-style vector collection" begin
    onebody = zeros(Float64, 4, 4)
    twobody = zeros(Float64, 4, 4, 4, 4)
    model = TwoBandModel(2, 2, 1, onebody, twobody)

    elist, mocoeffs = DiagonalizeAllMomentum(model, 1; numer_of_vectors=1, save=true)
    results = SolveAllSectors(model, SolverConfig(1; return_vectors=1))

    @test size(elist, 2) == 2
    @test Set(keys(mocoeffs)) == Set(result.momentum for result in results)
end

@testset "Unsupported spinful 6-index pipeline errors clearly" begin
    onebody = zeros(Float64, 4, 4)
    two6 = zeros(Float64, 2, 2, 2, 2, 2, 2)
    model = SpinfulMomentumModel(1, 1, 2, 2, onebody, two6)
    hilbertspace = BuildSector(model, 0)
    @test_throws ArgumentError BuildOperator(model, hilbertspace)
end
