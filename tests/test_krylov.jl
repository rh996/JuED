using KrylovKit
using LinearAlgebra

#linear operator

function f(v::Vector{Float64};x = 1)

    y = similar(v)
    for i in eachindex(v)
        y[i] = i*v[i]*x
    end
    return y
end




dim = 10
initv = rand(Float64, dim)


function g(v::Vector{Float64})
    
    return f(v, x=1)
end



# @show typeof(g)
lambda, u = eigsolve(g, dim,3; maxiter=1000, tol=1e-6,ishermitian=true)

@show lambda
@show u