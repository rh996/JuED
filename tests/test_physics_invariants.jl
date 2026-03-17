using LinearAlgebra
using Test

if isdefined(Main, :JuED)
    const EDMod = Main.JuED.EDMod
else
    include("../src/EDMain.jl")
    const EDMod = Main.EDMod
end

function spinless_particle_count(state::Integer, norbital::Int)
    count = 0
    for orbital in 1:norbital
        site = EDMod.FermionOperatorMod.basis_site_index(norbital, orbital)
        count += Int(((unsigned(state) >> (site - 1)) & 0x1) == 0x1)
    end
    return count
end

function spinless_2d_momentum(state::Integer, Nkx::Int, Nky::Int)
    norbital = Nkx * Nky
    momentum = 0
    for orbital in 1:norbital
        site = EDMod.FermionOperatorMod.basis_site_index(norbital, orbital)
        if ((unsigned(state) >> (site - 1)) & 0x1) == 0x1
            momentum = EDMod.MomentumHilbertSpace2DMod.momentum_add_2d(momentum, orbital - 1, Nkx, Nky)
        end
    end
    return momentum
end

@testset "2D spinless sector states preserve particle number and momentum" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    model = EDMod.SpinlessListModel(nparticle, Nkx, Nky, zeros(Float64, 4, 4), zeros(Float64, 4, 4, 4, 4))

    for momentum in 0:(Nkx * Nky - 1)
        hilbertspace = EDMod.BuildSector(model, momentum)
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

    list_model = EDMod.SpinlessListModel(2, 2, 2, list_onebody, list_twobody)
    list_hilbert = EDMod.BuildSector(list_model, 1)
    list_operator, _ = EDMod.BuildOperator(list_model, list_hilbert)
    @test isapprox(Matrix(list_operator), Matrix(list_operator)'; atol=1e-12)

    momentum_onebody = zeros(Float64, 4, 4)
    momentum_twobody = zeros(Float64, 2, 2, 2, 2, 2, 2)
    momentum_twobody[1, 1, 1, 1, 1, 1] = 0.5

    momentum_model = EDMod.SpinlessMomentumModel(2, 2, 2, momentum_onebody, momentum_twobody)
    momentum_hilbert = EDMod.BuildSector(momentum_model, 1)
    momentum_operator, _ = EDMod.BuildOperator(momentum_model, momentum_hilbert)
    @test isapprox(Matrix(momentum_operator), Matrix(momentum_operator)'; atol=1e-12)
end
