using SparseArrays




m=4
n = 4
nonzeros_val :: Vector{Int32} = [1,3,2,4,5,6,7,8]
nonzeros_row :: Vector{Int32} = [1,4,3,2,1,2,1,2]
nonzeros_col :: Vector{Int32} = [1,1,1,2,3,3,4,4]
indptr :: Vector{Int32}= [1,4,5,7,9]



A = SparseMatrixCSC{Int32,Int32}(m,n,indptr, nonzeros_row, nonzeros_val)




let
    temp = 1
    counting = 0
    ptr = [1]

    for i in [1,1,1,2,3,3,4,4]
        
        if temp ==i
            counting += 1
        else
            # counting +=1
            push!(ptr, counting+ptr[end])
            counting=1
            
        end

        temp = i

    end
    push!(ptr, counting+ptr[end])
    println(ptr)
end

# A = sparse(nonzeros_row, indptr, nonzeros_val, m, n)