@testset "2D spinless sector states preserve particle number and momentum" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    model = SpinlessListModel(nparticle, Nkx, Nky, zeros(Float64, 4, 4), zeros(Float64, 4, 4, 4, 4))

    for momentum in 0:(Nkx * Nky - 1)
        hilbertspace = BuildSector(model, momentum)
        for state in hilbertspace.hilbert
            @test spinless_particle_count(state, Nkx * Nky) == nparticle
            @test spinless_2d_momentum(state, Nkx, Nky) == momentum
        end
    end
end

@testset "Small representative Hamiltonians are Hermitian" begin
    list_onebody = Diagonal([1.0, -0.5, 0.75, 0.25]) |> Matrix
    list_twobody = zeros(Float64, 4, 4, 4, 4)
    list_twobody[1, 2, 2, 1] = 0.3
    list_twobody[2, 1, 1, 2] = 0.3
    list_twobody[3, 4, 4, 3] = -0.2
    list_twobody[4, 3, 3, 4] = -0.2

    list_model = SpinlessListModel(2, 2, 2, list_onebody, list_twobody)
    list_hilbert = BuildSector(list_model, 1)
    list_operator, _ = BuildOperator(list_model, list_hilbert)
    @test isapprox(Matrix(list_operator), Matrix(list_operator)'; atol=1e-12)

    momentum_onebody = zeros(Float64, 4, 4)
    momentum_twobody = zeros(Float64, 2, 2, 2, 2, 2, 2)
    momentum_twobody[1, 1, 1, 1, 1, 1] = 0.5

    momentum_model = SpinlessMomentumModel(2, 2, 2, momentum_onebody, momentum_twobody)
    momentum_hilbert = BuildSector(momentum_model, 1)
    momentum_operator, _ = BuildOperator(momentum_model, momentum_hilbert)
    @test isapprox(Matrix(momentum_operator), Matrix(momentum_operator)'; atol=1e-12)
end
