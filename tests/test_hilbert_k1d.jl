include("../src/EDMain.jl")
using .EDMod


hibertspace = EDMod.MomentumHilbertSpace1DMod.MomentumHilbertSpace1D{Int64}(1, 6, 1, [])

hilbert = EDMod.MomentumHilbertSpace1DMod.BuildHilbert(hibertspace)


println(length(hilbert))
println(hilbert)
for i in eachindex(hilbert)
    println(bitstring(hilbert[i])[end-5:end])
end
