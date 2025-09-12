include("../src/EDMain.jl")

using .EDMod

using HDF5
using LinearAlgebra

function test(nparticle,Nkx,Nky)
    ## Read the data
    filename = "./data/twoband_TMD_degree_3.89_gate_15.0_nkx_$(Nkx)_nky_$(Nky).h5"
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
    e,v = DiagonalizeOneMomentum(modelparams,0,5; returnvectors=2)
    @show size(e)
    @show size(v)

    @show typeof(v[1])
    

    rdm1 = EDMod.RDM1(modelparams,v[1],0)

    @show tr(rdm1)


    # @showtime EDMod.RDM2_cache(modelparams,v[1],0)
    # cache = "data/rdm2_cache_pnumber_6_K6x3_knumber_0.jld2"
    # cache = "./data/rdm2_cache_pnumber_8_K6x4_knumber_0.jld2"
    # cache = "./data/rdm2_chahe_K6x3.jld2"
    # @showtime rdm2 = EDMod.RDM2(modelparams,v[1],0,cache)

    @showtime rdm2 = EDMod.RDM2(modelparams,v[1],0)

    sum = 0.0 + 0.0im
    for i in 1:Nkx*Nky
        for j in 1:Nkx*Nky
            sum+= rdm2[i,j,i,j]
        end
    end

    @show sum

    comparison = zeros(ComplexF64,Nkx*Nky,Nkx*Nky)
    for i in 1:Nkx*Nky
        for j in 1:Nkx*Nky
            for k in 1:Nkx*Nky
                comparison[i,j] += rdm2[i,k,j,k]/(nparticle-1)

            end
    
        end
    end

    @show norm(comparison-rdm1)


    @showtime rdm3 = EDMod.RDM3(modelparams,v[1],1)
    # @showtime rdm3 = EDMod.RDM3_single_thread(modelparams,v[1],1)

    size(rdm3)

end

test(8,6,4)