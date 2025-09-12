include("../src/SpinMomentumHilbertSpace1D.jl")
using .SpinMomentumHilbertSpace1DMod

nalpha = 2
nbeta = 2
orbital = 4
k = 0
hilbertspace = SpinMomentumHilbertSpace1DMod.SpinMomentumHilbertSpace1D{Int64}(nalpha,nbeta,orbital,k,[])
hilbert = SpinMomentumHilbertSpace1DMod.BuildSpinHilbert(hilbertspace)


println(length(hilbert))
println(hilbert)
for i in eachindex(hilbert)
    println(bitstring(hilbert[i])[end-orbital*2+1:end])
end
