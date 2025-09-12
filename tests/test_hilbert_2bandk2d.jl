include("../src/EDMain.jl")

using .EDMod


nparticle = 6
Nkx = 3
Nky = 6

let 
    result = 0
    for k in 0:(Nkx*Nky-1)
        println("k = ", k)
        hibertspace = EDMod.TwoBandMomentumHilbertSpace2DMod.TwoBandMomentumHilbertSpace2D{Int64}(nparticle, Nkx, Nky, k, [])
        hilbert =  EDMod.TwoBandMomentumHilbertSpace2DMod.BuildTwoBandHilbert(hibertspace)
        curr = length(hilbert)
        result += curr
        println("length = ", curr)
    end
    println("total length = ", result)
end
