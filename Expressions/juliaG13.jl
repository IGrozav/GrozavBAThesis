function G_13s(p::Number, z::Number)::Complex{Float64}
    return (-im)*sqrt(p)*sqrt(1 - z^2)*(M(p-mass^2/4-im*mass*sqrt(p)*z) - M(p-mass^2/4+im*mass*sqrt(p)*z))
end