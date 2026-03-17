module PublicAPIMod

export DiagonalizeOneMomentum, InputModel, InputTwoBandModel, DiagonalizeAllMomentum, HamiltonianAction
export SpinlessListModel, SpinlessMomentumModel, SpinfulListModel, SpinfulMomentumModel, TwoBandModel
export SolverConfig, BuildSector, BuildOperator, SolveSector, SolveAllSectors
export RDMWorkspace, CompactRDM2, CompactRDM3, RDM1, RDM2, RDM3, RDM2Compact, RDM3Compact, RDM2_cache, todense
export compact_rdm_filename, save_compact_rdm, load_compact_rdm2, load_compact_rdm3
export BasisSpaces

using SparseArrays
using KrylovKit

import ..InternalMod
using ..InternalMod.ModelTypesMod
using ..InternalMod.DensityMatricesMod
using ..InternalMod.HamiltonianConstructorMod: HamiltonianAction
using ..InternalMod.IndexTypesMod: choose_state_type

const BasisSpaces = InternalMod.BasisSpaces

struct SolverConfig
    neigen::Int
    return_vectors::Int
    matrixfree::Bool
    tol::Float64
    maxiter::Int
    which::Symbol
    ishermitian::Bool
end

function SolverConfig(
    neigen::Int;
    return_vectors::Int=1,
    matrixfree::Bool=false,
    tol::Real=1e-6,
    maxiter::Int=1000,
    which::Symbol=:SR,
    ishermitian::Bool=true,
)
    return SolverConfig(neigen, return_vectors, matrixfree, Float64(tol), maxiter, which, ishermitian)
end

SpinlessListModel(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T} =
    ModelParams2DSpinlessList{T}(nparticle, Nkx, Nky, OneBody, TwoBody)

SpinfulListModel(nalpha::Int, nbeta::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T} =
    ModelParams2DSpinList{T}(nalpha, nbeta, Nkx, Nky, OneBody, TwoBody)

SpinlessMomentumModel(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,6}) where {T} =
    ModelParams2DSpinless{T}(nparticle, Nkx, Nky, OneBody, TwoBody)

SpinfulMomentumModel(nalpha::Int, nbeta::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,6}) where {T} =
    ModelParams2DSpin{T}(nalpha, nbeta, Nkx, Nky, OneBody, TwoBody)

TwoBandModel(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T} =
    ModelParams2DTwoBand{T}(nparticle, Nkx, Nky, OneBody, TwoBody)

InputModel(args...) = _input_model(args...)
_input_model(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T} =
    SpinlessListModel(nparticle, Nkx, Nky, OneBody, TwoBody)
_input_model(nalpha::Int, nbeta::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T} =
    SpinfulListModel(nalpha, nbeta, Nkx, Nky, OneBody, TwoBody)
_input_model(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,6}) where {T} =
    SpinlessMomentumModel(nparticle, Nkx, Nky, OneBody, TwoBody)
_input_model(nalpha::Int, nbeta::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,6}) where {T} =
    SpinfulMomentumModel(nalpha, nbeta, Nkx, Nky, OneBody, TwoBody)

InputTwoBandModel(nparticle::Int, Nkx::Int, Nky::Int, OneBody::Array{T,2}, TwoBody::Array{T,4}) where {T} =
    TwoBandModel(nparticle, Nkx, Nky, OneBody, TwoBody)

@inline _sector_count(model) = model.Nkx * model.Nky
@inline _model_value_type(model) = eltype(model.OneBody)
@inline _operator_value_type(model) = eltype(model.OneBody)

@inline _state_bits(model::Union{ModelParams2DSpinlessList,ModelParams2DSpinless}) = model.Nkx * model.Nky
@inline _state_bits(model::Union{ModelParams2DSpinList,ModelParams2DTwoBand}) = model.Nkx * model.Nky * 2
@inline _state_bits(model::ModelParams2DSpin) = model.Nkx * model.Nky * 2

@inline function _sector_hilbertspace(model::Union{ModelParams2DSpinlessList,ModelParams2DSpinless}, momentum::Int, ::Type{Ti}) where {Ti}
    return BasisSpaces.MomentumHilbertSpace2D{Ti}(model.nparticle, model.Nkx, model.Nky, momentum, Ti[])
end

@inline function _sector_hilbertspace(model::Union{ModelParams2DSpinList,ModelParams2DSpin}, momentum::Int, ::Type{Ti}) where {Ti}
    return BasisSpaces.SpinMomentumHilbertSpace2D{Ti}(model.nalpha, model.nbeta, model.Nkx, model.Nky, momentum, Ti[])
end

@inline function _sector_hilbertspace(model::ModelParams2DTwoBand, momentum::Int, ::Type{Ti}) where {Ti}
    return BasisSpaces.TwoBandMomentumHilbertSpace2D{Ti}(model.nparticle, model.Nkx, model.Nky, momentum, Ti[])
end

function BuildSector(model, momentum::Int; use_cache::Bool=true)
    indtype = choose_state_type(_state_bits(model))
    hilbertspace = _sector_hilbertspace(model, momentum, indtype)
    BasisSpaces.build_hilbert!(hilbertspace; use_cache)
    return hilbertspace
end

@inline _hamiltonian_constructor(model::ModelParams2DSpinless) = InternalMod.HamiltonianConstructorMod.Hamiltonian_momentum_constructor
@inline _hamiltonian_constructor(model::Union{ModelParams2DSpinlessList,ModelParams2DSpinList,ModelParams2DTwoBand}) = InternalMod.HamiltonianConstructorMod.Hamiltonian_list_constructor

function _sparse_operator_from_constructor(constructor, model, hilbertspace)
    data, row, indptr, dim = constructor(model.OneBody, model.TwoBody, hilbertspace)
    pointertype = eltype(row)
    operator = SparseMatrixCSC{_operator_value_type(model),pointertype}(dim, dim, indptr, row, data)
    return operator, dim
end

function BuildOperator(model::Union{ModelParams2DSpinless,ModelParams2DSpinlessList,ModelParams2DSpinList,ModelParams2DTwoBand}, hilbertspace; matrixfree::Bool=false)
    if matrixfree
        return InternalMod.HamiltonianConstructorMod.HamiltonianAction(model.OneBody, model.TwoBody, hilbertspace)
    end
    return _sparse_operator_from_constructor(_hamiltonian_constructor(model), model, hilbertspace)
end

function BuildOperator(model::ModelParams2DSpin, hilbertspace; matrixfree::Bool=false)
    throw(ArgumentError("Spinful 6-index models are not yet supported by the reusable ED solve pipeline."))
end

function _build_sector_problem(model, momentum::Int, config::SolverConfig; use_cache::Bool=true)
    hilbertspace = BuildSector(model, momentum; use_cache)
    operator, dim = BuildOperator(model, hilbertspace; matrixfree=config.matrixfree)
    return (; momentum, hilbertspace, operator, dim, matrixfree=config.matrixfree)
end

function _solve_operator(operator::SparseMatrixCSC, dim::Int, config::SolverConfig)
    return eigsolve(operator, dim, config.neigen, config.which; maxiter=config.maxiter, tol=config.tol, ishermitian=config.ishermitian)
end

function _solve_operator(action, dim::Int, config::SolverConfig)
    return eigsolve(action, dim, config.neigen, config.which; maxiter=config.maxiter, tol=config.tol, ishermitian=config.ishermitian)
end

function SolveSector(model, momentum::Int, config::SolverConfig=SolverConfig(1); use_cache::Bool=true)
    problem = _build_sector_problem(model, momentum, config; use_cache)
    dim = problem.dim
    nreq = min(config.neigen, dim)
    if nreq == 0
        empty_values = Vector{_model_value_type(model)}()
        vectors = config.return_vectors == 0 ? nothing : Any[]
        return (; momentum, dim, values=empty_values, vectors, info=nothing, hilbertspace=problem.hilbertspace)
    end
    vals, vecs, info = _solve_operator(problem.operator, dim, SolverConfig(nreq; return_vectors=config.return_vectors, matrixfree=config.matrixfree, tol=config.tol, maxiter=config.maxiter, which=config.which, ishermitian=config.ishermitian))
    vectors = config.return_vectors == 0 ? nothing : vecs[1:min(config.return_vectors, length(vecs))]
    return (; momentum, dim, values=vals[1:nreq], vectors, info, hilbertspace=problem.hilbertspace)
end

function SolveAllSectors(model, config::SolverConfig=SolverConfig(1); use_cache::Bool=true)
    results = Vector{Any}(undef, _sector_count(model))
    for (ind, momentum) in enumerate(0:(_sector_count(model) - 1))
        results[ind] = SolveSector(model, momentum, config; use_cache)
    end
    return results
end

@inline _nan_fill(::Type{T}) where {T<:Real} = T(NaN)
@inline _nan_fill(::Type{Complex{T}}) where {T<:Real} = Complex{T}(NaN, NaN)

function _stack_sector_values(results, neigen::Int)
    T = promote_type(map(result -> eltype(result.values), results)...)
    matrix = fill(_nan_fill(T), neigen, length(results))
    for (col, result) in enumerate(results)
        matrix[1:length(result.values), col] = result.values
    end
    return matrix
end

function _saved_vectors(results)
    mocoeffs = Dict{Int,Matrix}()
    for result in results
        result.vectors === nothing && continue
        mocoeffs[result.momentum] = hcat(result.vectors...)
    end
    return mocoeffs
end

function DiagonalizeOneMomentum(model::Union{ModelParams2DSpinless,ModelParams2DSpinlessList,ModelParams2DSpinList,ModelParams2DTwoBand}, momentum::Int, neigenv::Int; returnvectors::Int=1, matrixfree::Bool=false)
    result = SolveSector(model, momentum, SolverConfig(neigenv; return_vectors=returnvectors, matrixfree))
    if returnvectors == 0
        return result.values
    end
    return result.values, result.vectors
end

function _diagonalize_all_momentum(model, neigenv::Int; return_vectors::Int=0, save::Bool=false, matrixfree::Bool=false)
    results = SolveAllSectors(model, SolverConfig(neigenv; return_vectors, matrixfree))
    elist = _stack_sector_values(results, neigenv)
    if save
        return elist, _saved_vectors(results)
    end
    return elist
end

function DiagonalizeAllMomentum(model::Union{ModelParams2DSpinless,ModelParams2DSpinList}, neigenv::Int; matrixfree::Bool=false)
    return _diagonalize_all_momentum(model, neigenv; return_vectors=0, matrixfree)
end

function DiagonalizeAllMomentum(model::ModelParams2DSpinlessList, neigenv::Int; numer_of_vectors::Int=3, save::Bool=false, matrixfree::Bool=false)
    return_vectors = save ? numer_of_vectors : 0
    return _diagonalize_all_momentum(model, neigenv; return_vectors, save, matrixfree)
end

function DiagonalizeAllMomentum(model::ModelParams2DTwoBand, neigenv::Int; numer_of_vectors::Int=1, save::Bool=false, matrixfree::Bool=false)
    return_vectors = save ? numer_of_vectors : 0
    return _diagonalize_all_momentum(model, neigenv; return_vectors, save, matrixfree)
end

end
