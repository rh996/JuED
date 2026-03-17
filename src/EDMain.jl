module EDMod

include("Internal.jl")
include("PublicAPI.jl")

export DiagonalizeOneMomentum, InputModel, InputTwoBandModel, DiagonalizeAllMomentum, HamiltonianAction
export SpinlessListModel, SpinlessMomentumModel, SpinfulListModel, SpinfulMomentumModel, TwoBandModel
export SolverConfig, BuildSector, BuildOperator, SolveSector, SolveAllSectors
export RDMWorkspace, CompactRDM2, CompactRDM3, RDM1, RDM2, RDM3, RDM2Compact, RDM3Compact, RDM2_cache, todense
export compact_rdm_filename, save_compact_rdm, load_compact_rdm2, load_compact_rdm3
export BasisSpaces

using .PublicAPIMod
using .InternalMod.HamiltonianConstructorMod: HamiltonianAction, Hamiltonian_list_constructor, Hamiltonian_momentum_constructor
using .InternalMod.DensityMatricesMod: RDM3_single_thread, RDM2_naive, RDM3_naive

const BasisSpaces = PublicAPIMod.BasisSpaces

const ModelTypesMod = InternalMod.ModelTypesMod
const MomentumUtilsMod = InternalMod.MomentumUtilsMod
const BasisBuildersMod = InternalMod.BasisBuildersMod
const IndexTypesMod = InternalMod.IndexTypesMod
const FermionOperatorMod = InternalMod.FermionOperatorMod
const HilbertSpaceMod = InternalMod.HilbertSpaceMod
const MomentumHilbertSpace1DMod = InternalMod.MomentumHilbertSpace1DMod
const MomentumHilbertSpace2DMod = InternalMod.MomentumHilbertSpace2DMod
const SpinHilbertSpaceMod = InternalMod.SpinHilbertSpaceMod
const SpinMomentumHilbertSpace1DMod = InternalMod.SpinMomentumHilbertSpace1DMod
const SpinMomentumHilbertSpace2DMod = InternalMod.SpinMomentumHilbertSpace2DMod
const TwoBandMomentumHilbertSpace2DMod = InternalMod.TwoBandMomentumHilbertSpace2DMod
const BasisSpacesMod = InternalMod.BasisSpacesMod
const HamiltonianConstructorMod = InternalMod.HamiltonianConstructorMod
const DensityMatricesMod = InternalMod.DensityMatricesMod

end
