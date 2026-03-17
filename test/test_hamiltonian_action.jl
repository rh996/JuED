@testset "List Hamiltonian action matches sparse operator" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    momentum = 1

    onebody = Diagonal([1.0, 2.0, 3.0, 4.0]) |> Matrix
    twobody = zeros(Float64, 4, 4, 4, 4)
    model = SpinlessListModel(nparticle, Nkx, Nky, onebody, twobody)

    hilbertspace = BuildSector(model, momentum)
    sparse_h, dim = BuildOperator(model, hilbertspace)
    action, action_dim = BuildOperator(model, hilbertspace; matrixfree=true)

    @test dim == action_dim

    v = collect(1.0:dim)
    @test sparse_h * v == action(v)

    sparse_result = SolveSector(model, momentum, SolverConfig(1; return_vectors=0, matrixfree=false))
    action_result = SolveSector(model, momentum, SolverConfig(1; return_vectors=0, matrixfree=true))
    @test isapprox(sparse_result.values[1], action_result.values[1]; atol=1e-8)
end

@testset "Momentum Hamiltonian action matches sparse operator" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    momentum = 1

    onebody = zeros(Float64, 4, 4)
    twobody = zeros(Float64, Nkx, Nky, Nkx, Nky, Nkx, Nky)
    twobody[1, 1, 1, 1, 1, 1] = 0.5
    model = SpinlessMomentumModel(nparticle, Nkx, Nky, onebody, twobody)

    hilbertspace = BuildSector(model, momentum)
    sparse_h, dim = BuildOperator(model, hilbertspace)
    action, action_dim = BuildOperator(model, hilbertspace; matrixfree=true)

    @test dim == action_dim

    v = collect(1.0:dim)
    @test sparse_h * v == action(v)

    sparse_result = SolveSector(model, momentum, SolverConfig(1; return_vectors=0, matrixfree=false))
    action_result = SolveSector(model, momentum, SolverConfig(1; return_vectors=0, matrixfree=true))
    @test isapprox(sparse_result.values[1], action_result.values[1]; atol=1e-8)
end
