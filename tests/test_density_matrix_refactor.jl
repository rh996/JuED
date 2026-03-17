using LinearAlgebra
using Random
using Test
using JLD2

include("../src/EDMain.jl")
using .EDMod

const TEST_RDM3_LEFT_PERMS = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
const TEST_RDM3_RIGHT_PERMS = ((4, 5, 6), (4, 6, 5), (5, 4, 6), (5, 6, 4), (6, 4, 5), (6, 5, 4))

function find_nonempty_workspace(model)
    nsectors = model.Nkx * model.Nky
    for momentum in 0:(nsectors - 1)
        workspace = EDMod.RDMWorkspace(model, momentum)
        if !isempty(workspace.hilbert)
            return workspace
        end
    end
    error("No nonempty momentum sector found.")
end

function random_state(dim::Int)
    coeffs = randn(ComplexF64, dim)
    coeffs ./= norm(coeffs)
    return coeffs
end

function antisymmetrize_rdm2_local(rdm2)
    result = copy(rdm2)
    result .-= permutedims(result, (2, 1, 3, 4))
    result .-= permutedims(result, (1, 2, 4, 3))
    return result
end

function permutation_sign_local(perm)
    sign = 1
    for i in eachindex(perm)
        for j in (i + 1):length(perm)
            if perm[i] > perm[j]
                sign *= -1
            end
        end
    end
    return sign
end

function antisymmetrize_rdm3_local(rdm3)
    result = zeros(eltype(rdm3), size(rdm3)...)
    for eta in TEST_RDM3_LEFT_PERMS
        for sigma in TEST_RDM3_RIGHT_PERMS
            sign = permutation_sign_local(eta) * permutation_sign_local(sigma)
            result .+= sign .* permutedims(rdm3, (eta..., sigma...))
        end
    end
    return result
end

@testset "Spinless RDM workspace matches optimized and naive paths" begin
    Random.seed!(1234)
    onebody = zeros(ComplexF64, 4, 4)
    twobody = zeros(ComplexF64, 4, 4, 4, 4)
    model = EDMod.InputModel(2, 2, 2, onebody, twobody)
    workspace = find_nonempty_workspace(model)
    coeffs = random_state(length(workspace.hilbert))

    rdm1_workspace = EDMod.RDM1(workspace, coeffs)
    rdm1_model = EDMod.RDM1(model, coeffs, workspace.momentum)
    rdm2_workspace = EDMod.RDM2(workspace, coeffs)
    rdm2_compact = EDMod.RDM2Compact(workspace, coeffs)
    rdm2_compact_public = EDMod.RDM2(workspace, coeffs; representation=:compact)
    rdm2_naive = EDMod.RDM2_naive(workspace, coeffs)

    @test isapprox(rdm1_workspace, rdm1_model; atol=1e-10)
    @test isapprox(rdm2_workspace, EDMod.todense(rdm2_compact); atol=1e-10)
    @test isapprox(EDMod.todense(rdm2_compact_public), rdm2_workspace; atol=1e-10)
    @test length(rdm2_compact.elements) == length(workspace.rdm2_elements)
    for element in workspace.rdm2_elements
        i, j, k, l = element
        if (i, j) != (l, k)
            @test isapprox(rdm2_workspace[element...], rdm2_naive[element...]; atol=1e-10)
        end
    end

    cache_file = tempname() * ".jld2"
    try
        EDMod.RDM2_cache(workspace; file=cache_file)
        rdm2_compact_cached = EDMod.RDM2Compact(workspace, coeffs, cache_file)
        rdm2_cached = EDMod.RDM2(workspace, coeffs, cache_file)
        rdm2_cached_public = EDMod.RDM2(workspace, coeffs, cache_file; representation=:compact)
        @test isapprox(EDMod.todense(rdm2_compact_cached), rdm2_workspace; atol=1e-10)
        @test isapprox(rdm2_workspace, rdm2_cached; atol=1e-10)
        @test isapprox(EDMod.todense(rdm2_cached_public), rdm2_workspace; atol=1e-10)

        payload = load(cache_file)
        save(
            cache_file,
            "schema_version", payload["schema_version"],
            "cache_kind", payload["cache_kind"],
            "model_kind", payload["model_kind"],
            "Nkx", payload["Nkx"],
            "Nky", payload["Nky"],
            "momentum", Int32(workspace.momentum + 1),
            "nparticle", payload["nparticle"],
            "entries", payload["entries"],
            "k_$(workspace.momentum)", payload["entries"],
        )
        @test_throws ArgumentError EDMod.RDM2(workspace, coeffs, cache_file)
    finally
        isfile(cache_file) && rm(cache_file; force=true)
    end
end

@testset "Spinless RDM3 workspace matches single-thread and naive references" begin
    Random.seed!(5678)
    onebody = zeros(ComplexF64, 5, 5)
    twobody = zeros(ComplexF64, 5, 5, 5, 5)
    model = EDMod.InputModel(3, 5, 1, onebody, twobody)
    workspace = find_nonempty_workspace(model)
    coeffs = random_state(length(workspace.hilbert))

    rdm3_threaded = EDMod.RDM3(workspace, coeffs)
    rdm3_compact = EDMod.RDM3Compact(workspace, coeffs)
    rdm3_compact_public = EDMod.RDM3(workspace, coeffs; representation=:compact)
    rdm3_single = EDMod.RDM3_single_thread(workspace, coeffs)
    rdm3_naive = EDMod.RDM3_naive(workspace, coeffs)

    @test isapprox(rdm3_threaded, EDMod.todense(rdm3_compact); atol=1e-10)
    @test isapprox(EDMod.todense(rdm3_compact_public), rdm3_threaded; atol=1e-10)
    @test length(rdm3_compact.elements) == length(workspace.rdm3_elements)
    @test isapprox(rdm3_threaded, rdm3_single; atol=1e-10)
    for element in workspace.rdm3_elements
        i, j, k, l, m, n = element
        if (i, j, k) != (l, m, n)
            @test isapprox(rdm3_threaded[element...], rdm3_naive[element...]; atol=1e-10)
        end
    end
end

@testset "Two-band workspace preserves public RDM wrappers" begin
    Random.seed!(9012)
    onebody = zeros(ComplexF64, 4, 4)
    twobody = zeros(ComplexF64, 4, 4, 4, 4)
    model = EDMod.InputTwoBandModel(2, 2, 1, onebody, twobody)
    workspace = find_nonempty_workspace(model)
    coeffs = random_state(length(workspace.hilbert))

    @test isapprox(EDMod.RDM1(workspace, coeffs), EDMod.RDM1(model, coeffs, workspace.momentum); atol=1e-10)
    @test isapprox(EDMod.RDM2(workspace, coeffs), EDMod.RDM2(model, coeffs, workspace.momentum); atol=1e-10)
    @test isapprox(EDMod.todense(EDMod.RDM2(model, coeffs, workspace.momentum; representation=:compact)), EDMod.RDM2(workspace, coeffs); atol=1e-10)
    @test isapprox(EDMod.RDM2(workspace, coeffs), EDMod.todense(EDMod.RDM2Compact(workspace, coeffs)); atol=1e-10)
end
