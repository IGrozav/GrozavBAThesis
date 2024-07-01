function G_34s(p::Number, z::Number)::Complex{Float64}
    return (-((mass - (2*im)*sqrt(p)*z)*M(p-mass^2/4-im*mass*sqrt(p)*z)) + (mass + (2*im)*sqrt(p)*z)*M(p-mass^2/4+im*mass*sqrt(p)*z))/2
end