include("../src/EDMain.jl")

using .EDMod
using HDF5
using LinearAlgebra
using PythonCall
plt = pyimport("matplotlib.pyplot")

function test(nalpha,nbeta,Nkx,Nky)
    filename = "/Users/runhou/Library/Mobile Documents/com~apple~CloudDocs/Project_icloud/JuED/data/SpinfulInteraction_3band_degree_3.5_gate_20.0_nkx_$(Nkx)_nky_$(Nky).h5"
    
    interactions = h5read(filename, "SpinfulInteraction")
    interactions_dnup = h5read(filename, "SpinfulInteraction_dnup")
    interactions_updn = h5read(filename, "SpinfulInteraction_updn")
    interactions_dn = h5read(filename, "SpinfulInteraction_dn")
    band = h5read(filename, "band")


    eri = interactions[1,1,1,1,:,:,:,:,:,:]
    eri_dn = interactions_dn[1,1,1,1,:,:,:,:,:,:]
    eri_updn = interactions_updn[1,1,1,1,:,:,:,:,:,:]
    eri_dnup = interactions_dnup[1,1,1,1,:,:,:,:,:,:]

    Eri = zeros(ComplexF64,Nkx*Nky*2,Nkx*Nky*2,Nkx*Nky*2,Nkx*Nky*2)
    Onebody = zeros(ComplexF64,Nkx*Nky*2,Nkx*Nky*2)

    #onebody matrix
    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            ik = ikx + iky*Nkx

            #for another valley
            mikx = mod(-ikx,Nkx)
            miky = mod(-iky,Nky)

            Onebody[2*ik+1,2*ik+1] = band[ikx+1,iky+1,1]
            Onebody[2*ik+2,2*ik+2] = band[mikx+1,miky+1,1]
        end
    end
    
    #two body matrix
    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            for ikpx in 0:Nkx-1
                for ikpy in 0:Nky-1
                    for iqx in 0:Nkx-1
                        for iqy in 0:Nky-1
                            s1 = (mod(ikx + iqx,Nkx) + mod(iky + iqy,Nky)*Nkx)*2 
                            s2 = (ikx + iky*Nkx)*2
                            s3 = (mod(ikpx-iqx,Nkx) + mod(ikpy-iqy,Nky)*Nkx)*2 
                            s4 = (ikpx + ikpy*Nkx)*2 
                            Eri[s1+1,s3+1,s2+1,s4+1] += eri[ikx+1,iky+1,ikpx+1,ikpy+1,iqx+1,iqy+1]
                            Eri[s1+2,s3+2,s2+2,s4+2] += eri_dn[ikx+1,iky+1,ikpx+1,ikpy+1,iqx+1,iqy+1]
                            Eri[s1+1,s3+2,s2+1,s4+2] += eri_updn[ikx+1,iky+1,ikpx+1,ikpy+1,iqx+1,iqy+1]
                            Eri[s1+2,s3+1,s2+2,s4+1] += eri_dnup[ikx+1,iky+1,ikpx+1,ikpy+1,iqx+1,iqy+1]
                        end
                    end
                end
            end
        end
        
    end

    #check hermitian
    @show norm(Eri - permutedims(conj.(Eri), (3,4,1,2)))


    modelparams = InputModel(nalpha,nbeta,Nkx,Nky,Onebody,Eri)
    # modelparams2 = InputModel(nbeta,nalpha,Nkx,Nky,Onebody,Eri)
    elist = DiagonalizeAllMomentum(modelparams, 5)
    # elist2 = DiagonalizeAllMomentum(modelparams2, 5)

    # plot the energy

    X = [[i for i in 1:Nkx*Nky] for j in 1:5]
    X = hcat(X...)
    
    plt.scatter(X, elist', s=40, alpha=1.0 ,c="blue", marker="x")
    # plt.scatter(X, elist2', s=40, alpha=0.5, c="red", marker="o")
    plt.show()
    plt.close()


end


let 
    nalpha = 12
    nbeta = 0
    Nkx = 6
    Nky = 3
    test(nalpha,nbeta,Nkx,Nky)
    
end