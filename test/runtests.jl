using Test
using JuED

@testset "JuED" begin
    include("../tests/test_refactor_hilbert_spaces.jl")
    include("../tests/test_hilbert_kspin1d.jl")
    include("../tests/test_fermion_operator.jl")
    include("../tests/test_hamiltonian_action.jl")
    include("../tests/test_physics_invariants.jl")
    include("../tests/test_density_matrix.jl")
    include("../tests/test_density_matrix_refactor.jl")
    include("../tests/test_phase5_api.jl")
end
