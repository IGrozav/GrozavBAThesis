function G_11s(p::Number, z::Number)::Complex{Float64}
    return mass^2/4 + sqrt(p)^2 + M(p-mass^2/4-im*mass*sqrt(p)*z)*M(p-mass^2/4+im*mass*sqrt(p)*z)
end