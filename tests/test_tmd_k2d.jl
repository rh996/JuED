
include("../src/EDMain.jl")

using .EDMod
using HDF5

using PythonCall
plt = pyimport("matplotlib.pyplot")



function test(nparticle,Nkx,Nky)
    # filename = "./data/twoband_TMD_degree_3.89_gate_15.0_nkx_$(Nkx)_nky_$(Nky).h5"
    filename = "./data/twoband_TMD_degree_3.5_gate_10.0_nkx_$(Nkx)_nky_$(Nky).h5"
    interactions = h5read(filename, "interactions")
    band = h5read(filename, "band")

    Eri = interactions[1,1,1,1,:,:,:,:,:,:]
    Onebody = zeros(ComplexF64,Nkx*Nky,Nkx*Nky)
    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            ik = ikx + iky*Nkx
            Onebody[ik+1,ik+1] = band[ikx+1,iky+1,1]
        end
    end

 
    modelparams = InputModel(nparticle,Nkx,Nky,Onebody,Eri)
    elist = DiagonalizeAllMomentum(modelparams, 5)

 
    X = [[i for i in 1:Nkx*Nky] for j in 1:5]
    X = hcat(X...)
    
    plt.scatter(X, elist')
    plt.show()
    plt.close()

end



@time test(9,3,9)