include("../src/MomentumHilbertSpace1D.jl")
using .MomentumHilbertSpace1DMod


hibertspace = MomentumHilbertSpace1DMod.MomentumHilbertSpace1D{Int64}(1, 6, 1, [])

hilbert = MomentumHilbertSpace1DMod.BuildHilbert(hibertspace)


println(length(hilbert))
println(hilbert)
for i in eachindex(hilbert)
    println(bitstring(hilbert[i])[end-5:end])
end

