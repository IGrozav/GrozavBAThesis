function K_32ss(p::Number, q::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (-2*(sqrt(p)*z - sqrt(q)*zs)*(sqrt(p)*sqrt(1 - z^2) - sqrt(q)*y*sqrt(1 - zs^2)))/(sqrt(p)^2 + sqrt(q)^2 - 2*sqrt(p)*sqrt(q)*(z*zs + y*sqrt(1 - z^2)*sqrt(1 - zs^2)))
end