using Test

if isdefined(Main, :JuED)
    const EDMod = Main.JuED.EDMod
    const BasisSpaces = Main.JuED.BasisSpaces
else
    include("../src/EDMain.jl")
    const EDMod = Main.EDMod
    const BasisSpaces = Main.EDMod.BasisSpaces
end

@testset "Shared momentum utilities" begin
    @test BasisSpaces.momentum_add_1d(3, 4, 5) == 2
    @test BasisSpaces.momentum_sub_1d(1, 3, 5) == 3
    @test BasisSpaces.momentum_add_2d(1, 3, 2, 2) == 2
    @test BasisSpaces.momentum_sub_2d(0, 3, 2, 2) == 3
    @test BasisSpaces.momentum_add_2d(2, 3, 3, 2) == BasisSpaces.momentum_add_2d(2, 3, 3, 2)
end

@testset "2D momentum sectors partition the spinless Hilbert space" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    dims = Int[]
    for k in 0:(Nkx * Nky - 1)
        hilbertspace = BasisSpaces.MomentumHilbertSpace2D{Int32}(nparticle, Nkx, Nky, k, Int32[])
        push!(dims, length(BasisSpaces.build_hilbert!(hilbertspace)))
    end
    @test sum(dims) == binomial(Nkx * Nky, nparticle)
end

@testset "Momentum basis cache policy preserves exact contents" begin
    hilbertspace = BasisSpaces.MomentumHilbertSpace2D{Int32}(2, 2, 2, 1, Int32[])
    cached = BasisSpaces.build_hilbert!(hilbertspace; use_cache=true)
    uncached = BasisSpaces.build_hilbert!(hilbertspace; use_cache=false)
    @test cached == uncached
    @test cached == Int32[3, 12]
end

@testset "2D spin-momentum sectors partition the fixed-spin Hilbert space" begin
    nalpha = 1
    nbeta = 1
    Nkx = 2
    Nky = 2
    dims = Int[]
    for k in 0:(Nkx * Nky - 1)
        hilbertspace = BasisSpaces.SpinMomentumHilbertSpace2D{Int32}(nalpha, nbeta, Nkx, Nky, k, Int32[])
        push!(dims, length(BasisSpaces.build_hilbert!(hilbertspace)))
    end
    @test sum(dims) == binomial(Nkx * Nky, nalpha) * binomial(Nkx * Nky, nbeta)
end

@testset "Spin Hilbert builder respects cache policy" begin
    spin_hilbertspace = BasisSpaces.SpinHilbertSpace(1, 1, 2, Int64[])
    cached = BasisSpaces.build_hilbert!(spin_hilbertspace; use_cache=true)
    uncached = BasisSpaces.build_hilbert!(spin_hilbertspace; use_cache=false)
    @test cached == uncached
    @test cached == Int64[3, 9, 6, 12]
    @test spin_hilbertspace.hilbert == cached
end

@testset "2D two-band sectors partition the full Hilbert space" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    dims = Int[]
    for k in 0:(Nkx * Nky - 1)
        hilbertspace = BasisSpaces.TwoBandMomentumHilbertSpace2D{Int32}(nparticle, Nkx, Nky, k, Int32[])
        push!(dims, length(BasisSpaces.build_hilbert!(hilbertspace)))
    end
    @test sum(dims) == binomial(2 * Nkx * Nky, nparticle)
end
