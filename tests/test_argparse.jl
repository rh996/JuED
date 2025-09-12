include("../src/EDMain.jl")
using PythonCall
plt = pyimport("matplotlib.pyplot")
using HDF5

using .EDMod

function test(nparticle,Nkx,Nky)
    ## Read the data
    filename = "./data/twoband_TMD_degree_3.5_gate_10.0_nkx_$(Nkx)_nky_$(Nky).h5"
    interactions = h5read(filename, "interactions")
    band = h5read(filename, "band")

    Eri = interactions[1,1,1,1,:,:,:,:,:,:]
    Onebody = band[:,:,1]


    Eri_matrix = zeros(ComplexF64,Nkx*Nky,Nkx*Nky,Nkx*Nky,Nkx*Nky)
    Onebody_matrix = zeros(ComplexF64,Nkx*Nky,Nkx*Nky)
    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            for ikpx in 0:Nkx-1
                for ikpy in 0:Nky-1
                    for iqx in 0:Nkx-1
                        for iqy in 0:Nky-1
                            s1 = mod(ikx + iqx,Nkx) + mod(iky + iqy,Nky)*Nkx
                            s2 = ikx + iky*Nkx
                            s3 = mod(ikpx-iqx,Nkx) + mod(ikpy-iqy,Nky)*Nkx
                            s4 = ikpx + ikpy*Nkx
                            Eri_matrix[s1+1,s3+1,s2+1,s4+1] = Eri[ikx+1,iky+1,ikpx+1,ikpy+1,iqx+1,iqy+1]

                        end
                    end
                end
            end
        end
    end


    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            ik = ikx + iky*Nkx
            Onebody_matrix[ik+1,ik+1] = Onebody[ikx+1,iky+1]
        end
    end



    # test the diagonalization
    modelparams = InputModel(nparticle,Nkx,Nky,Onebody_matrix,Eri_matrix)
    elist,vector = DiagonalizeAllMomentum(modelparams, 5;save=true)

    # for k in 0:Nkx*Nky-1
    #     h5write("test.h5","vec_k_$k",vector[k])
    # end




    # e,v = DiagonalizeOneMomentum(modelparams,1,5;returnvectors=3)

    # @show e

    # v_array = hcat(v...)
    # @show typeof(v_array)
    # @show size(v_array)


    # h5write("test.h5", "v", v_array)
    # h5write("test.h5", "e", e)



    # plot the energy

    X = [[i for i in 1:Nkx*Nky] for j in 1:5]
    X = hcat(X...)
    
    plt.scatter(X, elist')
    plt.show()
    plt.close()

end

if length(ARGS) != 3
    println("Usage: julia test_argparse.jl particle Nkx Nky")
    exit(1)
end
particle  = parse(Int,ARGS[1])
Nkx = parse(Int,ARGS[2])
Nky = parse(Int,ARGS[3])

println("particle: ",particle)
println("Nkx: ",Nkx)
println("Nky: ",Nky)
test(particle,Nkx,Nky)