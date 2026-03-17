using Test

include("../src/FermionOperator.jl")
using .FermionOperatorMod

@testset "Basis site mapping documents reverse orbital order" begin
    @test basis_site_index(6, 1) == 6
    @test basis_site_index(6, 6) == 1
    @test basis_site_index(4, 2) == 3
end

@testset "Pure fermion kernels preserve integer type and sign conventions" begin
    creation = creation_kernel(Int32(0b0010), 1)
    @test creation == (Int32(0b0011), Int8(-1))
    @test typeof(first(creation)) == Int32

    annihilation = annihilation_kernel(Int32(0b1010), 2)
    @test annihilation == (Int32(0b1000), Int8(-1))
    @test typeof(first(annihilation)) == Int32

    @test creation_kernel(Int32(0b0010), 2) === nothing
    @test annihilation_kernel(Int32(0b0010), 1) === nothing

    string_result = apply_operator_string(Int32(0b0110), (1,), (2,))
    @test string_result == (Int32(0b0101), Int8(1))
end

@testset "Mutable wrappers match the pure kernels" begin
    fermion = FermionOperator(Int32(0b0110), 1)
    AnnihilationOperator!(fermion, 2)
    CreationOperator!(fermion, 1)

    @test fermion.state == Int32(0b0101)
    @test fermion.fermion_sign == Int8(1)

    failed = FermionOperator(Int64(0b0001), 1)
    CreationOperator!(failed, 1)
    @test failed.state == Int64(0b0001)
    @test failed.fermion_sign == Int8(0)
end

@testset "Signed state widths reject unsupported top-bit sites" begin
    high32 = creation_kernel(Int32(0), 31)
    @test high32 == (Int32(1) << 30, Int8(1))
    @test_throws ArgumentError creation_kernel(Int32(0), 32)

    high64 = creation_kernel(Int64(0), 63)
    @test high64 == (Int64(1) << 62, Int8(1))
    @test_throws ArgumentError creation_kernel(Int64(0), 64)
end
