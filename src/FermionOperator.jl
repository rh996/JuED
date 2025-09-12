
module FermionOperatorMod

export FermionOperator, CreationOperator!, AnnihilationOperator!

mutable struct FermionOperator
    state :: Int64
    fermion_sign :: Int64
end


function CreationOperator!(fermion_operator::FermionOperator, isite::Int64)
    # if isite < 0 || isite > 63
    #     throw(ArgumentError("isite must be between 0 and 63"))
    # end

    count = fermion_operator.state
    count = count >> (isite - 1)
    is_odd = 0

    if count & 1 != 0
        fermion_operator.fermion_sign = 0
    else
        fermion_operator.state = fermion_operator.state | (1 << (isite - 1))
        while count!= 0
            count &= count - 1
            is_odd ⊻= 1
        end
        if is_odd==1
            fermion_operator.fermion_sign *= -1
        else
            fermion_operator.fermion_sign *= 1
        end
    end
end

function AnnihilationOperator!(fermion_operator::FermionOperator, isite::Int64)
    # if isite < 0 || isite > 63
    #     throw(ArgumentError("isite must be between 0 and 63"))
    # end

    count = fermion_operator.state
    count = count >> (isite - 1)
    is_odd = 0

    if count & 1 != 1
        fermion_operator.fermion_sign = 0
    else
        fermion_operator.state = fermion_operator.state & (~(1 << (isite - 1)))
        while count!= 0
            count &= count - 1
            is_odd ⊻= 1
        end
        if is_odd==0
            fermion_operator.fermion_sign *= -1
        else
            fermion_operator.fermion_sign *= 1
        end
    end
    
end
end


