include("../src/HilbertSpace.jl")
# include("../src/FermionOperator.jl")
# include("../src/HamiltonianConstructor.jl")
include("../src/EDMain.jl")
using .EDMod

using Arpack
using SparseArrays
# using .HilbertSpaceMod

#two site hubbard module
hilbertspace = HilbertSpaceMod.GeneralHilbertSpace{Int128}(2, 6, [])
hibertarr = HilbertSpaceMod.BuildHilbert(hilbertspace)

#hopping matrix
t = -1.0
hop = zeros(Float64, 6, 6)

spinhop = [0.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 0.0]

hop[1:2:end, 1:2:end] = t * spinhop
hop[2:2:end, 2:2:end] = t * spinhop

# @show hop

#onsite interaction
U = 4.0
onsite = zeros(Float64, 6, 6, 6, 6)
onsite[1, 2, 1, 2] = U / 2
onsite[2, 1, 2, 1] = U / 2
onsite[3, 4, 3, 4] = U / 2
onsite[4, 3, 4, 3] = U / 2
onsite[5, 6, 5, 6] = U / 2
onsite[6, 5, 6, 5] = U / 2

test = EDMod.HamiltonianConstructorMod.GeneralHamiltonian(hop, onsite, hilbertspace)


# @show Matrix(test)
# data

e, lambda = eigs(test, nev=14, which=:SR)

println(e)

Matrix(test)

# using LinearAlgebra
# hamiltonian = Matrix(test)

# e2,l2 = eigen(hamiltonian)

# print(e2)
# dense = Matrix(hamiltonian)
