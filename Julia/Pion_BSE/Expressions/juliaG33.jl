function G_33s(p::Number, z::Number)::Complex{Float64}
    return -1/4*mass^2 + sqrt(p)^2*(1 - 2*z^2) + M(p-mass^2/4-im*mass*sqrt(p)*z)*M(p-mass^2/4+im*mass*sqrt(p)*z)
end