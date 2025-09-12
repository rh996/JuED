module TwoBandMomentumHilbertSpace2DMod

using ..HilbertSpaceMod

mutable struct TwoBandMomentumHilbertSpace2D{Ti}<:AbstractHilbertSpace
    nparticle::Int64
    Nkx::Int64
    Nky::Int64
    momentum::Int64
    hilbert::Array{Ti,1}
end


@inline function momentum_add_2d(k1,k2,Nkx,Nky)
    k1x = mod(k1,Nkx)
    k1y = fld(k1,Nkx)
    k2x = mod(k2,Nkx)
    k2y = fld(k2,Nkx)
    return mod((k1x+k2x),Nkx) + mod((k1y+k2y),Nky)*Nkx
    
end

@inline function momentum_sub_2d(k1,k2,Nkx,Nky)
    k1x = mod(k1,Nkx)
    k1y = fld(k1,Nkx)
    k2x = mod(k2,Nkx)
    k2y = fld(k2,Nkx)
    return mod((k1x-k2x),Nkx) + mod((k1y-k2y),Nky)*Nkx
    
end


function _dfs(ne::Int64,no::Int64,k_curr::Int64,hilbertspace::TwoBandMomentumHilbertSpace2D)
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    systemsize = Nkx*Nky

    if ne>no
        return []
    end
    if ne==no
        k_new = 0
        for i in 0:no-1
            k1 = (systemsize - 1-i)
            k_new = momentum_add_2d(k_new,k1,Nkx,Nky)
        end
        if k_new == k_curr
            a::Int64 = 0
            for i in 1:ne
                a |= 1<<((i-1)*2)
            end
            return [a]
        else
            return []
        end
    elseif ne==0
        if k_curr == 0
            return [0]
        else
            return []
        end
    else
        kind = systemsize-no
        k1 = k_curr
        k2 = kind
        k_new = momentum_sub_2d(k1,k2,Nkx,Nky)
        
        left = _dfs(ne,no-1,k_curr,hilbertspace)
        right = _dfs(ne-1,no-1,k_new,hilbertspace)
        shifted_right = right .+ (1<<((no-1)*2))
        curr = vcat(left,shifted_right)
        return curr
    end
end

function BuildHilbert(nparticle,momentum,hilbertspace::TwoBandMomentumHilbertSpace2D)
    norbital = hilbertspace.Nkx*hilbertspace.Nky
    k = momentum
    hilbert = _dfs(nparticle,norbital,k,hilbertspace)
    return hilbert
    
end

function BuildTwoBandHilbert(hilbertspace::TwoBandMomentumHilbertSpace2D)
    nparticle = hilbertspace.nparticle
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    norbital = hilbertspace.Nkx*hilbertspace.Nky
    momentum = hilbertspace.momentum
    result :: typeof(hilbertspace.hilbert)= []

    for nalpha in 0:nparticle
        nbeta = nparticle - nalpha
        for k1 in 0:norbital-1
            for k2 in 0:norbital-1
                k_tot = momentum_add_2d(k1,k2,Nkx,Nky)
                if k_tot == momentum
                    hilbert_alpha = BuildHilbert(nalpha,k1,hilbertspace)
                    hilbert_beta = BuildHilbert(nbeta,k2,hilbertspace)
                    for i in eachindex(hilbert_beta)
                        hilbert_beta[i] = hilbert_beta[i] << 1
                    end
                    hilbert_tot = zeros(typeof(hilbertspace.hilbert).parameters[1],length(hilbert_alpha)*length(hilbert_beta))
                    for i in eachindex(hilbert_alpha)
                        for j in eachindex(hilbert_beta)
                            hilbert_tot[(i-1)*length(hilbert_beta)+j] = hilbert_alpha[i] | hilbert_beta[j]
                        end
                    end
                    result = vcat(result,hilbert_tot)
                end
            end
        end
    end




    hilbertspace.hilbert = result
    return result
end





end