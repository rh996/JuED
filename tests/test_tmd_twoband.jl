include("../src/EDMain.jl")

using .EDMod
using HDF5
using LinearAlgebra
using PythonCall

plt = pyimport("matplotlib.pyplot")
function test(nparticle,Nkx,Nky)
    filename = "./data/twoband_TMD_degree_3.5_gate_20.0_nkx_$(Nkx)_nky_$(Nky).h5"
    interactions = h5read(filename, "interactions")
    band = h5read(filename, "band")

    nband = 2
    @show nband

    Eri_matrix = zeros(ComplexF64,Nkx*Nky*nband,Nkx*Nky*nband,Nkx*Nky*nband,Nkx*Nky*nband)
    Onebody_matrix = zeros(ComplexF64,Nkx*Nky*nband,Nkx*Nky*nband)

    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            for ikpx in 0:Nkx-1
                for ikpy in 0:Nky-1
                    for iqx in 0:Nkx-1
                        for iqy in 0:Nky-1
                            s1 = (mod(ikx + iqx,Nkx) + mod(iky + iqy,Nky)*Nkx)*nband 
                            s2 = (ikx + iky*Nkx)*nband
                            s3 = (mod(ikpx-iqx,Nkx) + mod(ikpy-iqy,Nky)*Nkx)*nband 
                            s4 = (ikpx + ikpy*Nkx)*nband 
                            for iband in 1:nband, jband in 1:nband, kband in 1:nband, lband in 1:nband
                                Eri_matrix[s1+iband,s3+kband,s2+jband,s4+lband] += interactions[iband,jband,kband,lband,ikx+1,iky+1,ikpx+1,ikpy+1,iqx+1,iqy+1]
                            end
                        end
                    end
                end
            end
        end
    end

    Eri_matrix = Eri_matrix + permutedims(conj.(Eri_matrix), (3,4,1,2))
    Eri_matrix  -= permutedims(Eri_matrix, (1,2,4,3))
    Eri_matrix  = Eri_matrix./4




    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            ik = ikx + iky*Nkx
            for iband in 1:nband
                Onebody_matrix[ik*nband+iband,ik*nband+iband] = band[ikx+1,iky+1,iband]
            end
        end
    end

    # display(Onebody_matrix)


    TwoBandModel = EDMod.InputTwoBandModel(nparticle,Nkx,Nky,Onebody_matrix,Eri_matrix)
    elist,mocoeffs = EDMod.DiagonalizeAllMomentum(TwoBandModel, 5; numer_of_vectors=5, save=true)
    @show elist'
    X = [[i for i in 0:Nkx*Nky-1] for j in 1:5]
    X = hcat(X...)

    plt.scatter(X, elist')
    plt.show()
    plt.close()

    rdm2 = EDMod.RDM2(TwoBandModel, mocoeffs[4][:,1], 4)
    rdm2 += EDMod.RDM2(TwoBandModel, mocoeffs[8][:,1], 8)
    rdm2 += EDMod.RDM2(TwoBandModel, mocoeffs[0][:,1], 0)
    rdm2 = rdm2./3



    # @show rdm2[1,2,3,4]
    # @show rdm2[4,3,2,1]
    # Save rdm2 as numpy .npy file
    np = pyimport("numpy")
    np.save("rdm2_2band.npy", rdm2)
    # rdm1 = EDMod.RDM1(TwoBandModel, mocoeffs[4][:,1], 4)
    # rdm1 += EDMod.RDM1(TwoBandModel, mocoeffs[8][:,1], 8)
    # rdm1 = rdm1./2
    # @show diag(rdm1[1:2:end,1:2:end])
    # @show diag(rdm1[2:2:end,2:2:end])
    # @show tr(rdm1)
    # @show tr(rdm1[1:2:end,1:2:end])
    # @show tr(rdm1[2:2:end,2:2:end])

    # norbital = nband*Nkx*Nky
    # sum = 0.0 + 0.0im
    # for i in 1:norbital
    #     for j in 1:norbital
    #         sum+= rdm2[i,j,i,j]
    #     end
    # end

    # @show sum

    # comparison = zeros(ComplexF64,Nkx*Nky*2,Nkx*Nky*2)
    # for i in 1:norbital
    #     for j in 1:norbital
    #         for k in 1:norbital 
    #             comparison[i,j] += rdm2[i,k,j,k]/(nparticle-1)

    #         end
    
    #     end
    # end

    # @show norm(comparison-rdm1)


    # @show real.(diag(rdm1))
end

test(4,2,6)