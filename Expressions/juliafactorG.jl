function factorG(p::Number, z::Number)::Complex{Float64}
    return 1/(A(p-mass^2/4-im*mass*sqrt(p)*z)*A(p-mass^2/4+im*mass*sqrt(p)*z)*(-1/4*mass^2 + sqrt(p)^2 - im*mass*sqrt(p)*z + M(p-mass^2/4-im*mass*sqrt(p)*z)^2)*(-1/4*mass^2 + sqrt(p)^2 + im*mass*sqrt(p)*z + M(p-mass^2/4+im*mass*sqrt(p)*z)^2))
end