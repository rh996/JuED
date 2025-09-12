include("../src/EDMain.jl")

using .EDMod

hibertspace = EDMod.MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{Int64}(4, 2,6, 3, [])

hilbert =  EDMod.MomentumHilbertSpace2DMod.BuildHilbert(hibertspace)


println(length(hilbert))
# println(hilbert)
# for i in eachindex(hilbert)
#     println(bitstring(hilbert[i]))
# end

###########################################

