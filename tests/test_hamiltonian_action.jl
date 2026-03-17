using Test
using LinearAlgebra
using SparseArrays

include("../src/EDMain.jl")
using .EDMod

@testset "List Hamiltonian action matches sparse CSC" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    momentum = 1

    onebody = Diagonal([1.0, 2.0, 3.0, 4.0]) |> Matrix
    twobody = zeros(Float64, 4, 4, 4, 4)
    model = EDMod.InputModel(nparticle, Nkx, Nky, onebody, twobody)

    hilbertspace = EDMod.MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{Int32}(nparticle, Nkx, Nky, momentum, [])
    EDMod.MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)

    data, row, indptr, dim = EDMod.Hamiltonian_list_constructor(onebody, twobody, hilbertspace)
    action, action_dim = EDMod.HamiltonianAction(onebody, twobody, hilbertspace)

    @test dim == action_dim

    pointertype = eltype(row)
    sparse_h = SparseMatrixCSC{Float64,pointertype}(dim, dim, indptr, row, data)
    v = collect(1.0:dim)

    @test sparse_h * v == action(v)

    sparse_vals, _ = EDMod.DiagonalizeOneMomentum(model, momentum, 1)
    action_vals, _ = EDMod.DiagonalizeOneMomentum(model, momentum, 1; matrixfree=true)
    @test isapprox(sparse_vals[1], action_vals[1]; atol=1e-8)
end

@testset "Momentum Hamiltonian action matches sparse CSC" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    momentum = 1

    onebody = zeros(Float64, 4, 4)
    twobody = zeros(Float64, Nkx, Nky, Nkx, Nky, Nkx, Nky)
    twobody[1, 1, 1, 1, 1, 1] = 0.5
    model = EDMod.InputModel(nparticle, Nkx, Nky, onebody, twobody)

    hilbertspace = EDMod.MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{Int32}(nparticle, Nkx, Nky, momentum, [])
    EDMod.MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)

    data, row, indptr, dim = EDMod.Hamiltonian_momentum_constructor(onebody, twobody, hilbertspace)
    action, action_dim = EDMod.HamiltonianAction(onebody, twobody, hilbertspace)

    @test dim == action_dim

    pointertype = eltype(row)
    sparse_h = SparseMatrixCSC{Float64,pointertype}(dim, dim, indptr, row, data)
    v = collect(1.0:dim)

    @test sparse_h * v == action(v)

    sparse_vals, _ = EDMod.DiagonalizeOneMomentum(model, momentum, 1)
    action_vals, _ = EDMod.DiagonalizeOneMomentum(model, momentum, 1; matrixfree=true)
    @test isapprox(sparse_vals[1], action_vals[1]; atol=1e-8)
end
