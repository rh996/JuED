module ModelTypesMod

export ModelParams2DSpinlessList, ModelParams2DSpinless, ModelParams2DSpinList, ModelParams2DSpin, ModelParams2DTwoBand

struct ModelParams2DSpinlessList{T}
    nparticle::Int
    Nkx::Int
    Nky::Int
    OneBody::Array{T,2}
    TwoBody::Array{T,4}
end

struct ModelParams2DSpinless{T}
    nparticle::Int
    Nkx::Int
    Nky::Int
    OneBody::Array{T,2}
    TwoBody::Array{T,6}
end


struct ModelParams2DSpinList{T}
    nalpha::Int
    nbeta::Int
    Nkx::Int
    Nky::Int
    OneBody::Array{T,2}
    TwoBody::Array{T,4}
end

struct ModelParams2DSpin{T}
    nalpha::Int
    nbeta::Int
    Nkx::Int
    Nky::Int
    OneBody::Array{T,2}
    TwoBody::Array{T,6}
end

struct ModelParams2DTwoBand{T}
    nparticle::Int
    Nkx::Int
    Nky::Int
    OneBody::Array{T,2}
    TwoBody::Array{T,4}
end



end