module BasisSpacesMod

export AbstractHilbertSpace
export GeneralHilbertSpace, MomentumHilbertSpace1D, MomentumHilbertSpace2D
export SpinHilbertSpace, SpinMomentumHilbertSpace1D, SpinMomentumHilbertSpace2D
export TwoBandMomentumHilbertSpace2D
export build_hilbert!, state_index_map
export momentum_add_1d, momentum_sub_1d, momentum_add_2d, momentum_sub_2d

import ..HilbertSpaceMod
import ..MomentumHilbertSpace1DMod
import ..MomentumHilbertSpace2DMod
import ..SpinHilbertSpaceMod
import ..SpinMomentumHilbertSpace1DMod
import ..SpinMomentumHilbertSpace2DMod
import ..TwoBandMomentumHilbertSpace2DMod

using ..HilbertSpaceMod: AbstractHilbertSpace, GeneralHilbertSpace
using ..MomentumHilbertSpace1DMod: MomentumHilbertSpace1D, momentum_add_1d, momentum_sub_1d
using ..MomentumHilbertSpace2DMod: MomentumHilbertSpace2D, momentum_add_2d, momentum_sub_2d
using ..SpinHilbertSpaceMod: SpinHilbertSpace
using ..SpinMomentumHilbertSpace1DMod: SpinMomentumHilbertSpace1D
using ..SpinMomentumHilbertSpace2DMod: SpinMomentumHilbertSpace2D
using ..TwoBandMomentumHilbertSpace2DMod: TwoBandMomentumHilbertSpace2D

build_hilbert!(hilbertspace::GeneralHilbertSpace; use_cache::Bool=true) =
    HilbertSpaceMod.BuildHilbert(hilbertspace; use_cache)

build_hilbert!(hilbertspace::MomentumHilbertSpace1D; use_cache::Bool=true) =
    MomentumHilbertSpace1DMod.BuildHilbert(hilbertspace; use_cache)

build_hilbert!(hilbertspace::MomentumHilbertSpace2D; use_cache::Bool=true) =
    MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace; use_cache)

build_hilbert!(hilbertspace::SpinHilbertSpace; use_cache::Bool=true) =
    SpinHilbertSpaceMod.BuildSpinHilbert(hilbertspace; use_cache)

build_hilbert!(hilbertspace::SpinMomentumHilbertSpace1D; use_cache::Bool=true) =
    SpinMomentumHilbertSpace1DMod.BuildSpinHilbert(hilbertspace; use_cache)

build_hilbert!(hilbertspace::SpinMomentumHilbertSpace2D; use_cache::Bool=true) =
    SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace; use_cache)

build_hilbert!(hilbertspace::TwoBandMomentumHilbertSpace2D; use_cache::Bool=true) =
    TwoBandMomentumHilbertSpace2DMod.BuildTwoBandHilbert(hilbertspace; use_cache)

state_index_map(hilbertspace::AbstractHilbertSpace, ::Type{Tv}=Int) where {Tv<:Integer} =
    HilbertSpaceMod.ToDict(hilbertspace, Tv)

end
