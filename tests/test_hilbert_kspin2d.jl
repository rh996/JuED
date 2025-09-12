include("../src/EDMain.jl")
using .EDMod

nalpha = 2
nbeta = 2
Nkx = 4
Nky = 4
k = 0
orbital = Nkx * Nky


hilbertspace = EDMod.SpinMomentumHilbertSpace2DMod.SpinMomentumHilbertSpace2D{Int64}(nalpha, nbeta, Nkx, Nky, k, [])
hilbert = EDMod.SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)


println(length(hilbert))
# println(hilbert)
# for i in eachindex(hilbert)
#     println(bitstring(hilbert[i])[end-orbital*2+1:end])
# end
