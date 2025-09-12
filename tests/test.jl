import Tullio
import Base
using Arpack

# println(count_ones(5))

# println(bitstring(5))
# # println("\n")
# println(bitstring(10))

# println(bitstring(5 | 10))

include("../src/HilbertSpace.jl")
using .HilbertSpaceMod

hilbertspace = HilbertSpaceMod.GeneralHilbertSpace{Int64}(3, 5, [])



@time hibert0 = HilbertSpaceMod.BuildHilbert(hilbertspace)

println(hibert0)

println(length(hibert0))


for i in eachindex(hibert0)
    println(bitstring(hibert0[i]))
end



include("../src/SpinHilbertSpace.jl")
using .SpinHilbertSpaceMod
spin_hilbertspace = SpinHilbertSpaceMod.SpinHilbertSpace(3, 1, 5, [])

@time hilbert = SpinHilbertSpaceMod.BuildSpinHilbert(spin_hilbertspace)
# nalpha = spin_hilbertspace.nalpha
# nbeta = spin_hilbertspace.nbeta
# hilbert = BuildHilbert(nbeta,spin_hilbertspace)

println(length(hilbert))
println(hilbert)
for i in eachindex(hilbert)
    println(bitstring(hilbert[i]))
end

ind_dict = ToDict(hilbertspace)

include("../src/FermionOperator.jl")
using .FermionOperatorMod
println("test fermion operator")

fermion_operator = FermionOperator(0b101010, 1)
println(bitstring(fermion_operator.state))
CreationOperator!(fermion_operator, 1)

println(bitstring(fermion_operator.state))
println(fermion_operator.fermion_sign)

AnnihilationOperator!(fermion_operator, 4)

println(bitstring(fermion_operator.state))
println(fermion_operator.fermion_sign)

CreationOperator!(fermion_operator, 5)

println(bitstring(fermion_operator.state))
println(fermion_operator.fermion_sign)
