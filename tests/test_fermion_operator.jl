using Test

if isdefined(Main, :JuED)
    const FermionOperatorMod = Main.JuED.EDMod.FermionOperatorMod
else
    include("../src/FermionOperator.jl")
    const FermionOperatorMod = Main.FermionOperatorMod
end

@testset "Basis site mapping documents reverse orbital order" begin
    @test FermionOperatorMod.basis_site_index(6, 1) == 6
    @test FermionOperatorMod.basis_site_index(6, 6) == 1
    @test FermionOperatorMod.basis_site_index(4, 2) == 3
end

@testset "Pure fermion kernels preserve integer type and sign conventions" begin
    creation = FermionOperatorMod.creation_kernel(Int32(0b0010), 1)
    @test creation == (Int32(0b0011), Int8(-1))
    @test typeof(first(creation)) == Int32

    annihilation = FermionOperatorMod.annihilation_kernel(Int32(0b1010), 2)
    @test annihilation == (Int32(0b1000), Int8(-1))
    @test typeof(first(annihilation)) == Int32

    @test FermionOperatorMod.creation_kernel(Int32(0b0010), 2) === nothing
    @test FermionOperatorMod.annihilation_kernel(Int32(0b0010), 1) === nothing

    string_result = FermionOperatorMod.apply_operator_string(Int32(0b0110), (1,), (2,))
    @test string_result == (Int32(0b0101), Int8(1))
end

@testset "Mutable wrappers match the pure kernels" begin
    fermion = FermionOperatorMod.FermionOperator(Int32(0b0110), 1)
    FermionOperatorMod.AnnihilationOperator!(fermion, 2)
    FermionOperatorMod.CreationOperator!(fermion, 1)

    @test fermion.state == Int32(0b0101)
    @test fermion.fermion_sign == Int8(1)

    failed = FermionOperatorMod.FermionOperator(Int64(0b0001), 1)
    FermionOperatorMod.CreationOperator!(failed, 1)
    @test failed.state == Int64(0b0001)
    @test failed.fermion_sign == Int8(0)
end

@testset "Signed state widths reject unsupported top-bit sites" begin
    high32 = FermionOperatorMod.creation_kernel(Int32(0), 31)
    @test high32 == (Int32(1) << 30, Int8(1))
    @test_throws ArgumentError FermionOperatorMod.creation_kernel(Int32(0), 32)

    high64 = FermionOperatorMod.creation_kernel(Int64(0), 63)
    @test high64 == (Int64(1) << 62, Int8(1))
    @test_throws ArgumentError FermionOperatorMod.creation_kernel(Int64(0), 64)
end
