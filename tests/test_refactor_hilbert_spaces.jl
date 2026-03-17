using Test

if isdefined(Main, :JuED)
    const EDMod = Main.JuED.EDMod
else
    include("../src/EDMain.jl")
    const EDMod = Main.EDMod
end

@testset "Shared momentum utilities" begin
    @test EDMod.MomentumHilbertSpace1DMod.momentum_add_1d(3, 4, 5) == 2
    @test EDMod.MomentumHilbertSpace1DMod.momentum_sub_1d(1, 3, 5) == 3
    @test EDMod.MomentumHilbertSpace2DMod.momentum_add_2d(1, 3, 2, 2) == 2
    @test EDMod.MomentumHilbertSpace2DMod.momentum_sub_2d(0, 3, 2, 2) == 3
    @test EDMod.SpinMomentumHilbertSpace2DMod.momentum_add_2d(2, 3, 3, 2) ==
          EDMod.MomentumHilbertSpace2DMod.momentum_add_2d(2, 3, 3, 2)
end

@testset "2D momentum sectors partition the spinless Hilbert space" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    dims = Int[]
    for k in 0:(Nkx * Nky - 1)
        hilbertspace = EDMod.MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{Int32}(nparticle, Nkx, Nky, k, [])
        push!(dims, length(EDMod.MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)))
    end
    @test sum(dims) == binomial(Nkx * Nky, nparticle)
end

@testset "Momentum basis cache policy preserves exact contents" begin
    hilbertspace = EDMod.MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{Int32}(2, 2, 2, 1, [])
    cached = EDMod.MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace; use_cache=true)
    uncached = EDMod.MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace; use_cache=false)
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
        hilbertspace = EDMod.SpinMomentumHilbertSpace2DMod.SpinMomentumHilbertSpace2D{Int32}(nalpha, nbeta, Nkx, Nky, k, [])
        push!(dims, length(EDMod.SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)))
    end
    @test sum(dims) == binomial(Nkx * Nky, nalpha) * binomial(Nkx * Nky, nbeta)
end

@testset "Spin Hilbert builder respects cache policy" begin
    spin_hilbertspace = EDMod.SpinHilbertSpaceMod.SpinHilbertSpace(1, 1, 2, [])
    cached = EDMod.SpinHilbertSpaceMod.BuildSpinHilbert(spin_hilbertspace; use_cache=true)
    uncached = EDMod.SpinHilbertSpaceMod.BuildSpinHilbert(spin_hilbertspace; use_cache=false)
    @test cached == uncached
    @test cached == Int64[3, 9, 6, 12]
end

@testset "2D two-band sectors partition the full Hilbert space" begin
    nparticle = 2
    Nkx = 2
    Nky = 2
    dims = Int[]
    for k in 0:(Nkx * Nky - 1)
        hilbertspace = EDMod.TwoBandMomentumHilbertSpace2DMod.TwoBandMomentumHilbertSpace2D{Int32}(nparticle, Nkx, Nky, k, [])
        push!(dims, length(EDMod.TwoBandMomentumHilbertSpace2DMod.BuildTwoBandHilbert(hilbertspace)))
    end
    @test sum(dims) == binomial(2 * Nkx * Nky, nparticle)
end
