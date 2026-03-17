module InternalMod

include("ModelTypes.jl")
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
include("BasisSpaces.jl")
include("HamiltonianConstructor.jl")
include("DensityMatrices.jl")

const BasisSpaces = BasisSpacesMod

end
