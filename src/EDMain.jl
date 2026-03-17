


module EDMod

include("ModelTypes.jl")
using .ModelTypesMod

include("MomentumUtils.jl")
include("BasisBuilders.jl")
include("IndexTypes.jl")
include("FermionOperator.jl")
include("HilbertSpace.jl")
include("SpinMomentumHilbertSpace2D.jl")
include("SpinMomentumHilbertSpace1D.jl")
include("MomentumHilbertSpace2D.jl")
include("MomentumHilbertSpace1D.jl")
include("SpinHilbertSpace.jl")
include("TwoBandMomentumHilbertSpace2D.jl")
include("HamiltonianConstructor.jl")
include("DensityMatrices.jl")

export DiagonalizeOneMomentum, InputModel, DiagonalizeAllMomentum

using .HamiltonianConstructorMod
using .DensityMatricesMod

using SparseArrays
using KrylovKit
using .MomentumHilbertSpace2DMod
using .MomentumHilbertSpace1DMod
using .SpinMomentumHilbertSpace2DMod
using .SpinMomentumHilbertSpace1DMod
using .SpinHilbertSpaceMod
using .TwoBandMomentumHilbertSpace2DMod
using .IndexTypesMod: choose_state_type


function InputModel(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T}
    return ModelParams2DSpinlessList{T}(nparticle, Nkx, Nky, OneBody, TwoBody)
end


function InputModel(nalpha::Int, nbeta::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T}
    return ModelParams2DSpinList{T}(nalpha, nbeta, Nkx, Nky, OneBody, TwoBody)
end

function InputModel(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,6}) where {T}
    return ModelParams2DSpinless{T}(nparticle, Nkx, Nky, OneBody, TwoBody)
end

function InputModel(nalpha::Int, nbeta::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,6}) where {T}
    return ModelParams2DSpin{T}(nalpha, nbeta, Nkx, Nky, OneBody, TwoBody)
end
function InputTwoBandModel(nparticle, Nkx, Nky, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T}
    return ModelParams2DTwoBand{T}(nparticle, Nkx, Nky, OneBody, TwoBody)
end


function DiagonalizeOneMomentum(ModelParams2DSpinless::ModelParams2DSpinless{T}, momentum::Int, neigenv::Int; returnvectors::Int=1) where {T}

    nparticle = ModelParams2DSpinless.nparticle
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    OneBody = ModelParams2DSpinless.OneBody
    TwoBody = ModelParams2DSpinless.TwoBody

    indtype = choose_state_type(Nkx * Nky)


    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)

    data, row, indptr, dim = Hamiltonian_momentum_constructor(OneBody, TwoBody, hilbertspace)

    @show dim

    if neigenv > dim
        neigenv = dim
    end

    pointertype = typeof(row[1])
    vals, vecs, info = eigsolve(SparseMatrixCSC{T,pointertype}(dim, dim, indptr, row, data), dim, neigenv, :SR; maxiter=1000, tol=1e-6, ishermitian=true)

    if returnvectors == 0
        return vals[1:neigenv]

    else
        return vals[1:neigenv], vecs[1:returnvectors]
    end


end


function DiagonalizeOneMomentum(ModelParams2DSpinless::ModelParams2DSpinlessList{T}, momentum::Int, neigenv::Int; returnvectors::Int=1) where {T}

    nparticle = ModelParams2DSpinless.nparticle
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    OneBody = ModelParams2DSpinless.OneBody
    TwoBody = ModelParams2DSpinless.TwoBody

    indtype = choose_state_type(Nkx * Nky)


    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)

    data, row, indptr, dim = Hamiltonian_list_constructor(OneBody, TwoBody, hilbertspace)

    @show dim

    if neigenv > dim
        neigenv = dim
    end

    pointertype = typeof(row[1])
    println("pointertype", pointertype)
    vals, vecs, info = eigsolve(SparseMatrixCSC{T,pointertype}(dim, dim, indptr, row, data), dim, neigenv, :SR; maxiter=1000, tol=1e-6, ishermitian=true)
    @show info
    if returnvectors == 0
        return vals[1:neigenv]

    else
        return vals, vecs[1:returnvectors]
    end


end


function DiagonalizeOneMomentum(ModelParams2DSpin::ModelParams2DSpinList{T}, momentum::Int, neigenv::Int; returnvectors::Int=1) where {T}

    nalpha = ModelParams2DSpin.nalpha
    nbeta = ModelParams2DSpin.nbeta
    Nkx = ModelParams2DSpin.Nkx
    Nky = ModelParams2DSpin.Nky
    OneBody = ModelParams2DSpin.OneBody
    TwoBody = ModelParams2DSpin.TwoBody

    indtype = choose_state_type(Nkx * Nky * 2)


    hilbertspace = SpinMomentumHilbertSpace2DMod.SpinMomentumHilbertSpace2D{indtype}(nalpha, nbeta, Nkx, Nky, momentum, [])
    SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)

    data, row, indptr, dim = Hamiltonian_list_constructor(OneBody, TwoBody, hilbertspace)

    @show dim

    if neigenv > dim
        neigenv = dim
    end

    pointertype = typeof(row[1])
    vals, vecs, info = eigsolve(SparseMatrixCSC{T,pointertype}(dim, dim, indptr, row, data), dim, neigenv, :SR; maxiter=1000, tol=1e-6, ishermitian=true)

    if returnvectors == 0
        return vals[1:neigenv]

    else
        return vals[1:neigenv], vecs[1:returnvectors]
    end


end





function DiagonalizeAllMomentum(ModelParams2DSpinless::ModelParams2DSpinless{T}, neigenv::Int) where {T}
    nparticle = ModelParams2DSpinless.nparticle
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    OneBody = ModelParams2DSpinless.OneBody
    TwoBody = ModelParams2DSpinless.TwoBody

    indtype = choose_state_type(Nkx * Nky)


    elist = []
    for k in 0:Nkx*Nky-1
        hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, k, [])
        hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)

        data, row, indptr, dim = Hamiltonian_momentum_constructor(OneBody, TwoBody, hilbertspace)

        @show dim

        if neigenv > dim
            neigenv = dim
        end
        # vals,vecs = eigs(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), nev=5, which=:SR)

        pointertype = typeof(row[1])
        vals, vecs, info = eigsolve(SparseMatrixCSC{T,pointertype}(dim, dim, indptr, row, data), dim, neigenv, :SR; maxiter=1000, tol=1e-6, ishermitian=true)
        # vals, vecs, info = eigsolve(SparseMatrixCSC{ComplexF64,Int32}(dim, dim, indptr, row, data),ones(ComplexF64,dim),10,:SR, Lanczos())
        push!(elist, vals[1:neigenv])

    end
    elist = hcat(elist...)


    return elist
end


function DiagonalizeAllMomentum(ModelParams2DSpinless::ModelParams2DSpinlessList{T}, neigenv::Int; numer_of_vectors=3, save=false) where {T}
    nparticle = ModelParams2DSpinless.nparticle
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    OneBody = ModelParams2DSpinless.OneBody
    TwoBody = ModelParams2DSpinless.TwoBody

    indtype = choose_state_type(Nkx * Nky)





    elist = []
    mocoeffs = Dict()
    for k in 0:Nkx*Nky-1
        hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, k, [])
        hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)

        data, row, indptr, dim = Hamiltonian_list_constructor(OneBody, TwoBody, hilbertspace)
        # println("Total memory consumption: ", (Base.summarysize(data)+Base.summarysize(row)+Base.summarysize(indptr))/1024^2, " MB")
        @show dim

        if neigenv > dim
            neigenv = dim
        end
        # vals,vecs = eigs(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), nev=5, which=:SR)
        pointertype = typeof(row[1])
        vals, vecs, info = eigsolve(SparseMatrixCSC{T,pointertype}(dim, dim, indptr, row, data), dim, neigenv, :SR; maxiter=1000, tol=1e-6, ishermitian=true)
        # vals, vecs, info = eigsolve(SparseMatrixCSC{ComplexF64,Int32}(dim, dim, indptr, row, data),ones(ComplexF64,dim),10,:SR, Lanczos())
        push!(elist, vals[1:neigenv])
        if save == true
            mocoeffs[k] = hcat(vecs[1:numer_of_vectors]...)
        end

    end
    elist = hcat(elist...)

    if save == true
        return elist, mocoeffs
    else
        return elist
    end
end


function DiagonalizeAllMomentum(ModelParams2DSpin::ModelParams2DSpinList{T}, neigenv::Int) where {T}
    nalpha = ModelParams2DSpin.nalpha
    nbeta = ModelParams2DSpin.nbeta
    Nkx = ModelParams2DSpin.Nkx
    Nky = ModelParams2DSpin.Nky
    OneBody = ModelParams2DSpin.OneBody
    TwoBody = ModelParams2DSpin.TwoBody

    indtype = choose_state_type(Nkx * Nky * 2)





    elist = []
    for k in 0:Nkx*Nky-1
        hilbertspace = SpinMomentumHilbertSpace2DMod.SpinMomentumHilbertSpace2D{indtype}(nalpha, nbeta, Nkx, Nky, k, [])
        SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)

        data, row, indptr, dim = Hamiltonian_list_constructor(OneBody, TwoBody, hilbertspace)

        @show dim

        if neigenv > dim
            neigenv = dim
        end
        # vals,vecs = eigs(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), nev=5, which=:SR)
        pointertype = typeof(row[1])
        vals, vecs, info = eigsolve(SparseMatrixCSC{T,pointertype}(dim, dim, indptr, row, data), dim, neigenv, :SR; maxiter=1000, tol=1e-6, ishermitian=true)
        # vals, vecs, info = eigsolve(SparseMatrixCSC{ComplexF64,Int32}(dim, dim, indptr, row, data),ones(ComplexF64,dim),10,:SR, Lanczos())
        push!(elist, vals[1:neigenv])

    end
    elist = hcat(elist...)


    return elist
end


function DiagonalizeAllMomentum(ModelParams2DTwoBand::ModelParams2DTwoBand{T}, neigenv::Int; numer_of_vectors=1, save=false) where {T}
    nparticle = ModelParams2DTwoBand.nparticle
    Nkx = ModelParams2DTwoBand.Nkx
    Nky = ModelParams2DTwoBand.Nky
    OneBody = ModelParams2DTwoBand.OneBody
    TwoBody = ModelParams2DTwoBand.TwoBody

    indtype = choose_state_type(Nkx * Nky * 2)
    elist = []
    mocoeffs = Dict()
    for k in 0:Nkx*Nky-1
        hilbertspace = TwoBandMomentumHilbertSpace2DMod.TwoBandMomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, k, [])
        TwoBandMomentumHilbertSpace2DMod.BuildTwoBandHilbert(hilbertspace)

        data, row, indptr, dim = Hamiltonian_list_constructor(OneBody, TwoBody, hilbertspace)

        @show dim

        if neigenv > dim
            neigenv = dim
        end
        # vals,vecs = eigs(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), nev=5, which=:SR)
        pointertype = typeof(row[1])
        vals, vecs, info = eigsolve(SparseMatrixCSC{T,pointertype}(dim, dim, indptr, row, data), dim, neigenv, :SR; maxiter=1000, tol=1e-6, ishermitian=true)
        # vals, vecs, info = eigsolve(SparseMatrixCSC{ComplexF64,Int32}(dim, dim, indptr, row, data),ones(ComplexF64,dim),10,:SR, Lanczos())
        push!(elist, vals[1:neigenv])
        if save == true
            mocoeffs[k] = hcat(vecs[1:numer_of_vectors]...)
        end

    end
    elist = hcat(elist...)

    if save == true
        return elist, mocoeffs
    else
        return elist
    end
end


end
