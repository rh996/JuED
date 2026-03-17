using JuED
using JLD2
using LinearAlgebra
using Random
using SparseArrays
using Test

const BasisSpaces = JuED.BasisSpaces
const Compat = JuED.EDMod
const FermionOperatorInternal = Compat.FermionOperatorMod

function find_nonempty_workspace(model)
    for momentum in 0:(model.Nkx * model.Nky - 1)
        workspace = RDMWorkspace(model, momentum)
        if !isempty(workspace.hilbert)
            return workspace
        end
    end
    error("No nonempty momentum sector found.")
end

random_state(dim::Int) = normalize(randn(ComplexF64, dim))

function contract_rdm2(rdm2, nparticle::Int)
    norbital = size(rdm2, 1)
    contracted = zeros(eltype(rdm2), norbital, norbital)
    for i in 1:norbital
        for j in 1:norbital
            for k in 1:norbital
                contracted[i, j] += rdm2[i, k, j, k] / nparticle
            end
        end
    end
    return contracted
end

function spinful_1d_counts(state::Integer, norbital::Int)
    nalpha = 0
    nbeta = 0
    for orbital in 0:(norbital - 1)
        nalpha += Int(((unsigned(state) >> (2 * orbital)) & 0x1) == 0x1)
        nbeta += Int(((unsigned(state) >> (2 * orbital + 1)) & 0x1) == 0x1)
    end
    return nalpha, nbeta
end

function spinful_1d_momentum(state::Integer, norbital::Int)
    momentum = 0
    for orbital in 0:(norbital - 1)
        if ((unsigned(state) >> (2 * orbital)) & 0x1) == 0x1
            momentum = BasisSpaces.momentum_add_1d(momentum, orbital, norbital)
        end
        if ((unsigned(state) >> (2 * orbital + 1)) & 0x1) == 0x1
            momentum = BasisSpaces.momentum_add_1d(momentum, orbital, norbital)
        end
    end
    return momentum
end

@inline basis_site_from_orbital(norbital::Int, orbital::Int) = norbital - orbital + 1

function spinless_particle_count(state::Integer, norbital::Int)
    count = 0
    for orbital in 1:norbital
        site = basis_site_from_orbital(norbital, orbital)
        count += Int(((unsigned(state) >> (site - 1)) & 0x1) == 0x1)
    end
    return count
end

function spinless_2d_momentum(state::Integer, Nkx::Int, Nky::Int)
    norbital = Nkx * Nky
    momentum = 0
    for orbital in 1:norbital
        site = basis_site_from_orbital(norbital, orbital)
        if ((unsigned(state) >> (site - 1)) & 0x1) == 0x1
            momentum = BasisSpaces.momentum_add_2d(momentum, orbital - 1, Nkx, Nky)
        end
    end
    return momentum
end
