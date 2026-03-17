include("TestSupport.jl")

@testset "JuED" begin
    include("test_phase1_namespace.jl")
    include("test_refactor_hilbert_spaces.jl")
    include("test_hilbert_kspin1d.jl")
    include("test_fermion_operator.jl")
    include("test_hamiltonian_action.jl")
    include("test_physics_invariants.jl")
    include("test_density_matrix.jl")
    include("test_density_matrix_refactor.jl")
    include("test_phase5_api.jl")
end
