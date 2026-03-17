using Test

if isdefined(Main, :JuED)
    const EDMod = Main.JuED.EDMod
else
    include("../src/EDMain.jl")
    const EDMod = Main.EDMod
end

function spinful_1d_counts(state::Integer, norbital::Int)
    nalpha = 0
    nbeta = 0
    for orbital in 0:(norbital - 1)
        nalpha += Int(((unsigned(state) >> (2 * orbital)) & 0x1) == 0x1)
        nbeta += Int(((unsigned(state) >> (2 * orbital + 1)) & 0x1) == 0x1)
    end
    return nalpha, nbeta
end

function spinful_1d_momentum(state::Integer, norbital::Int)
    momentum = 0
    for orbital in 0:(norbital - 1)
        if ((unsigned(state) >> (2 * orbital)) & 0x1) == 0x1
            momentum = EDMod.SpinMomentumHilbertSpace1DMod.momentum_add_1d(momentum, orbital, norbital)
        end
        if ((unsigned(state) >> (2 * orbital + 1)) & 0x1) == 0x1
            momentum = EDMod.SpinMomentumHilbertSpace1DMod.momentum_add_1d(momentum, orbital, norbital)
        end
    end
    return momentum
end

@testset "1D spin-momentum basis exact contents and quantum numbers" begin
    nalpha = 2
    nbeta = 2
    norbital = 4
    momentum = 0

    hilbertspace = EDMod.SpinMomentumHilbertSpace1DMod.SpinMomentumHilbertSpace1D{Int64}(nalpha, nbeta, norbital, momentum, Int64[])
    hilbert = EDMod.SpinMomentumHilbertSpace1DMod.BuildSpinHilbert(hilbertspace)

    @test hilbert == Int64[51, 45, 135, 120, 210, 204, 30, 180, 75, 225]

    for state in hilbert
        counts = spinful_1d_counts(state, norbital)
        @test counts == (nalpha, nbeta)
        @test spinful_1d_momentum(state, norbital) == momentum
    end
end

@testset "1D spin-momentum sectors partition the fixed-spin Hilbert space" begin
    nalpha = 2
    nbeta = 2
    norbital = 4
    dims = Int[]

    for momentum in 0:(norbital - 1)
        hilbertspace = EDMod.SpinMomentumHilbertSpace1DMod.SpinMomentumHilbertSpace1D{Int64}(nalpha, nbeta, norbital, momentum, Int64[])
        push!(dims, length(EDMod.SpinMomentumHilbertSpace1DMod.BuildSpinHilbert(hilbertspace)))
    end

    @test sum(dims) == binomial(norbital, nalpha) * binomial(norbital, nbeta)
end
