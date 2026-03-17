@testset "1D spin-momentum basis exact contents and quantum numbers" begin
    nalpha = 2
    nbeta = 2
    norbital = 4
    momentum = 0

    hilbertspace = BasisSpaces.SpinMomentumHilbertSpace1D{Int64}(nalpha, nbeta, norbital, momentum, Int64[])
    hilbert = BasisSpaces.build_hilbert!(hilbertspace)

    @test hilbert == Int64[51, 45, 135, 120, 210, 204, 30, 180, 75, 225]

    for state in hilbert
        @test spinful_1d_counts(state, norbital) == (nalpha, nbeta)
        @test spinful_1d_momentum(state, norbital) == momentum
    end
end

@testset "1D spin-momentum sectors partition the fixed-spin Hilbert space" begin
    nalpha = 2
    nbeta = 2
    norbital = 4
    dims = Int[]

    for momentum in 0:(norbital - 1)
        hilbertspace = BasisSpaces.SpinMomentumHilbertSpace1D{Int64}(nalpha, nbeta, norbital, momentum, Int64[])
        push!(dims, length(BasisSpaces.build_hilbert!(hilbertspace)))
    end

    @test sum(dims) == binomial(norbital, nalpha) * binomial(norbital, nbeta)
end
