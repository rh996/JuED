@testset "Spinless RDM workspace matches compact, cache, and model wrappers" begin
    Random.seed!(1234)
    onebody = zeros(ComplexF64, 4, 4)
    twobody = zeros(ComplexF64, 4, 4, 4, 4)
    model = SpinlessListModel(2, 2, 2, onebody, twobody)
    workspace = find_nonempty_workspace(model)
    coeffs = random_state(length(workspace.hilbert))

    rdm1_workspace = RDM1(workspace, coeffs)
    rdm1_model = RDM1(model, coeffs, workspace.momentum)
    rdm2_workspace = RDM2(workspace, coeffs)
    rdm2_compact = RDM2Compact(workspace, coeffs)
    rdm2_compact_public = RDM2(workspace, coeffs; representation=:compact)
    rdm2_naive = Compat.RDM2_naive(workspace, coeffs)

    @test isapprox(rdm1_workspace, rdm1_model; atol=1e-10)
    @test isapprox(rdm2_workspace, todense(rdm2_compact); atol=1e-10)
    @test isapprox(todense(rdm2_compact_public), rdm2_workspace; atol=1e-10)
    @test length(rdm2_compact.elements) == length(workspace.rdm2_elements)
    for element in workspace.rdm2_elements
        i, j, k, l = element
        if (i, j) != (l, k)
            @test isapprox(rdm2_workspace[element...], rdm2_naive[element...]; atol=1e-10)
        end
    end

    cache_file = tempname() * ".jld2"
    compact_file = tempname() * ".jld2"
    try
        RDM2_cache(workspace; file=cache_file)
        rdm2_compact_cached = RDM2Compact(workspace, coeffs, cache_file)
        rdm2_cached = RDM2(workspace, coeffs, cache_file)
        rdm2_cached_public = RDM2(workspace, coeffs, cache_file; representation=:compact)
        @test isapprox(todense(rdm2_compact_cached), rdm2_workspace; atol=1e-10)
        @test isapprox(rdm2_workspace, rdm2_cached; atol=1e-10)
        @test isapprox(todense(rdm2_cached_public), rdm2_workspace; atol=1e-10)

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
        @test_throws ArgumentError RDM2(workspace, coeffs, cache_file)

        rdm2_compact_saved = RDM2Compact(workspace, coeffs; savefile=compact_file)
        rdm2_compact_loaded = load_compact_rdm2(workspace; file=compact_file)
        rdm2_compact_loaded_public = RDM2Compact(workspace; file=compact_file)
        @test isapprox(todense(rdm2_compact_loaded), rdm2_workspace; atol=1e-10)
        @test rdm2_compact_loaded.elements == rdm2_compact_saved.elements
        @test rdm2_compact_loaded.values == rdm2_compact_saved.values
        @test rdm2_compact_loaded_public.values == rdm2_compact_saved.values

        compact_payload = load(compact_file)
        save(compact_file, merge(compact_payload, Dict("momentum" => Int32(workspace.momentum + 1))))
        @test_throws ArgumentError load_compact_rdm2(workspace; file=compact_file)
    finally
        isfile(cache_file) && rm(cache_file; force=true)
        isfile(compact_file) && rm(compact_file; force=true)
    end

    model_compact_file = tempname() * ".jld2"
    try
        RDM2Compact(model, coeffs, workspace.momentum; savefile=model_compact_file)
        @test isapprox(todense(load_compact_rdm2(workspace; file=model_compact_file)), rdm2_workspace; atol=1e-10)
    finally
        isfile(model_compact_file) && rm(model_compact_file; force=true)
    end
end

@testset "Spinless RDM3 workspace matches compact, saved, and internal references" begin
    Random.seed!(5678)
    onebody = zeros(ComplexF64, 5, 5)
    twobody = zeros(ComplexF64, 5, 5, 5, 5)
    model = SpinlessListModel(3, 5, 1, onebody, twobody)
    workspace = find_nonempty_workspace(model)
    coeffs = random_state(length(workspace.hilbert))

    rdm3_threaded = RDM3(workspace, coeffs)
    rdm3_compact = RDM3Compact(workspace, coeffs)
    rdm3_compact_public = RDM3(workspace, coeffs; representation=:compact)
    rdm3_single = Compat.RDM3_single_thread(workspace, coeffs)
    rdm3_naive = Compat.RDM3_naive(workspace, coeffs)

    @test isapprox(rdm3_threaded, todense(rdm3_compact); atol=1e-10)
    @test isapprox(todense(rdm3_compact_public), rdm3_threaded; atol=1e-10)
    @test length(rdm3_compact.elements) == length(workspace.rdm3_elements)
    @test isapprox(rdm3_threaded, rdm3_single; atol=1e-10)
    for element in workspace.rdm3_elements
        i, j, k, l, m, n = element
        if (i, j, k) != (l, m, n)
            @test isapprox(rdm3_threaded[element...], rdm3_naive[element...]; atol=1e-10)
        end
    end

    compact_file = tempname() * ".jld2"
    try
        rdm3_compact_saved = RDM3Compact(workspace, coeffs; savefile=compact_file)
        rdm3_compact_loaded = load_compact_rdm3(workspace; file=compact_file)
        rdm3_compact_loaded_public = RDM3Compact(workspace; file=compact_file)
        @test isapprox(todense(rdm3_compact_loaded), rdm3_threaded; atol=1e-10)
        @test rdm3_compact_loaded.elements == rdm3_compact_saved.elements
        @test rdm3_compact_loaded.values == rdm3_compact_saved.values
        @test rdm3_compact_loaded_public.values == rdm3_compact_saved.values

        compact_payload = load(compact_file)
        save(compact_file, merge(compact_payload, Dict("rdm_order" => Int32(2))))
        @test_throws ArgumentError load_compact_rdm3(workspace; file=compact_file)
    finally
        isfile(compact_file) && rm(compact_file; force=true)
    end

    model_compact_file = tempname() * ".jld2"
    try
        RDM3Compact(model, coeffs, workspace.momentum; savefile=model_compact_file)
        @test isapprox(todense(load_compact_rdm3(workspace; file=model_compact_file)), rdm3_threaded; atol=1e-10)
    finally
        isfile(model_compact_file) && rm(model_compact_file; force=true)
    end
end

@testset "Two-band workspace preserves public RDM wrappers" begin
    Random.seed!(9012)
    onebody = zeros(ComplexF64, 4, 4)
    twobody = zeros(ComplexF64, 4, 4, 4, 4)
    model = TwoBandModel(2, 2, 1, onebody, twobody)
    workspace = find_nonempty_workspace(model)
    coeffs = random_state(length(workspace.hilbert))

    @test isapprox(RDM1(workspace, coeffs), RDM1(model, coeffs, workspace.momentum); atol=1e-10)
    @test isapprox(RDM2(workspace, coeffs), RDM2(model, coeffs, workspace.momentum); atol=1e-10)
    @test isapprox(todense(RDM2(model, coeffs, workspace.momentum; representation=:compact)), RDM2(workspace, coeffs); atol=1e-10)
    @test isapprox(RDM2(workspace, coeffs), todense(RDM2Compact(workspace, coeffs)); atol=1e-10)
end
