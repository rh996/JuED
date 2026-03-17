@testset "Package exports BasisSpaces but not EDMod" begin
    @test :BasisSpaces in names(JuED)
    @test !(:EDMod in names(JuED))
    @test isdefined(JuED, :EDMod)
    @test isdefined(JuED.EDMod, :PublicAPIMod)
    @test isdefined(JuED.EDMod, :InternalMod)
end

@testset "BasisSpaces provides a unified basis namespace" begin
    general = BasisSpaces.GeneralHilbertSpace{Int32}(2, 4, Int32[])
    general_basis = BasisSpaces.build_hilbert!(general)

    @test general_basis == Int32[3, 5, 6, 9, 10, 12]
    @test BasisSpaces.state_index_map(general, Int32)[Int32(6)] == Int32(3)

    sector = BasisSpaces.MomentumHilbertSpace2D{Int32}(2, 2, 2, 1, Int32[])
    sector_basis = BasisSpaces.build_hilbert!(sector; use_cache=false)

    @test sector_basis == Int32[3, 12]
    @test BasisSpaces.state_index_map(sector, Int32)[Int32(12)] == Int32(2)
end
