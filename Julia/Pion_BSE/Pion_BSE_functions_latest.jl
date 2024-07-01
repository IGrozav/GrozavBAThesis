
#integration schemes
function Newtonroots(f::Function, e::Number, m::Number) #f..anonymous function, a...lower integration bound
    #b...Upper integration bound, e...Accuracy tolerance
    #initial guess
    root = Vector{Float64}(undef, m)
    for i = 1:m
        root[i] = cos(pi * (i - 0.25) / (m + 0.5))
    end
    h = 10^(-8)

    #Calculating the roots

    for j = 1:m
        xold = root[j]
        diffcon = 1.0
        n = 0
        while diffcon >= e && n < 10000
            global xnew = xold - f(xold) / (f(xold + h) - f(xold - h)) * 2 * h #Newton Method
            diffcon = abs(xold - xnew)
            xold = xnew
            n += 1

        end

        root[j] = xnew
    end
    return root
end

function legendre_poly(n::Int, x)
    # Initialize the first two Legendre polynomials
    P_n_minus_1 = 1.0
    P_n = x

    # Use the recurrence relation to compute higher-order Legendre polynomials
    for i = 2:n+1
        P_n_plus_1 = ((2 * i - 1) * x * P_n - (i - 1) * P_n_minus_1) / i
        P_n_minus_1 = P_n
        P_n = P_n_plus_1
    end
    return P_n, P_n_minus_1
end

function Gauss_Legendre(a::Number, b::Number, m::Int)

    Lpoly(x) = legendre_poly(m::Int, x)
    LP_n(x)=Lpoly(x)[2]
    LP_np1(x)=Lpoly(x)[1]
    e = 1.0 * 10^(-13)
    K = Newtonroots(LP_n::Function, e::Number, m::Number) #Calculating the Roots of the Legendre Polynomials in the interval[-1,1]
    K=reverse(K)
    w_i = 2 .* (1 .- K .^ 2) ./ ((m+1) .^ 2 .* LP_np1.(K) .^ 2)
    x_i=((b .- a) .* K .+ a .+ b) ./ 2

    return x_i, w_i .* (b .- a) ./ 2
end

function Gauss_Legendre1(L::Real, m::Int)

    Lpoly(x) = legendre_poly(m::Int, x)
    LP_n(x) = Lpoly(x)[2]
    LP_np1(x) = Lpoly(x)[1]
    e = 1.0 * 10^(-15)
    K = Newtonroots(LP_n::Function, e::Number, m::Number) #Calculating the Roots of the Legendre Polynomials in the interval[-1,1]
    K = reverse(K)
    w_i = 2 .* (1 .- K .^ 2) ./ ((m + 1) .^ 2 .* LP_np1.(K) .^ 2)
    x_i = (K .+ 1) ./ (1 .- K .+ 2 ./ L )

    return x_i, w_i .* (2 .+ 2 ./ L ) ./ (1 .- K .+ 2 ./ L ) .^ 2

end

function Gauss_Chebyshev2(a::Number, b::Number, m::Int)
    K = Vector{Float64}(undef, m)
    K = cos.((1:m) ./ (m .+ 1) .* pi)
    #Creating the Weights w_i
    global w_i = Vector{Any}(undef, m)
    global w_i = pi ./ (m .+ 1) .* sin.((1:m) ./ (m .+ 1) .* pi) .^ 2
    return ((b .- a) .* reverse(K) .+ a .+ b) ./ 2, (b .- a) ./ 2 .* w_i
end
#-------------------------------------------------------------------------------------
#This Section is For the quark propagator dressing functions

#Self implemented exponentially scaled modified bessel function of the first kind for n=1 and n=2
function bessel_map(n::Number, x::Number)::Number
    if real(x) < 0
        x = -x
    end
    if (n == 1)
        if real(x) >= 3.75
            return Complex(x) .^ (-0.5) .* (0.39894228 .- 0.03988024 .* (x ./ 3.75) .^ (-1) .- 0.00362018 .* (x ./ 3.75) .^ (-2) .+ 0.00163801 .* (x ./ 3.75) .^ (-3) .- 0.01031555 .* (x ./ 3.75) .^ (-4) .+ 0.02282967 .* (x ./ 3.75) .^ (-5) .- 0.02895312 .* (x ./ 3.75) .^ (-6) .+ 0.01787654 .* (x ./ 3.75) .^ (-7) .- 0.00420059 .* (x ./ 3.75) .^ (-8))
        else
            return exp.(-abs(real.(x))) .* x .* (0.5 .+ 0.87890594 .* (x ./ 3.75) .^ 2 .+ 0.51498869 .* (x ./ 3.75) .^ 4 .+ 0.15084934 .* (x ./ 3.75) .^ 6 .+ 0.02658733 .* (x ./ 3.75) .^ 8 .+ 0.00301532 .* (x ./ 3.75) .^ 10 .+ 0.00032411 .* (x ./ 3.75) .^ 12)
        end
    end
    if (n == 2)
        if real(x) >= 3.75
            return Complex(x) .^ (-0.5) .* (0.39894228 .+ 0.01328582 .* (x ./ 3.75) .^ (-1) .+ 0.00225319 .* (x ./ 3.75) .^ (-2) .- 0.00157565 .* (x ./ 3.75) .^ (-3) .+ 0.00916281 .* (x ./ 3.75) .^ (-4) .- 0.02057706 .* (x ./ 3.75) .^ (-5) .+ 0.02635537 .* (x ./ 3.75) .^ (-6) .- 0.01647633 .* (x ./ 3.75) .^ (-7) .+ 0.00392377 .* (x ./ 3.75) .^ (-8)) .- 2 .* Complex(x) .^ (-1.5) .* (0.39894228 .- 0.03988024 .* (x ./ 3.75) .^ (-1) .- 0.00362018 .* (x ./ 3.75) .^ (-2) .+ 0.00163801 .* (x ./ 3.75) .^ (-3) .- 0.01031555 .* (x ./ 3.75) .^ (-4) .+ 0.02282967 .* (x ./ 3.75) .^ (-5) .- 0.02895312 .* (x ./ 3.75) .^ (-6) .+ 0.01787654 .* (x ./ 3.75) .^ (-7) .- 0.00420059 .* (x ./ 3.75) .^ (-8))
        else
            return exp.(-abs(real.(x))) .* ((1 .+ 3.5156229 .* (x ./ 3.75) .^ 2 .+ 3.0899424 .* (x ./ 3.75) .^ 4 .+ 1.2067492 .* (x ./ 3.75) .^ 6 .+ 0.2659732 .* (x ./ 3.75) .^ 8 .+ 0.0360768 .* (x ./ 3.75) .^ 10 .+ 0.0045813 .* (x ./ 3.75) .^ 12) .- 2 .* (1 ./ 2 .+ 0.87890594 .* (x ./ 3.75) .^ 2 .+ 0.51498869 .* (x ./ 3.75) .^ 4 .+ 0.15084934 .* (x ./ 3.75) .^ 6 .+ 0.02658733 .* (x ./ 3.75) .^ 8 .+ 0.00301532 .* (x ./ 3.75) .^ 10 .+ 0.00032411 .* (x ./ 3.75) .^ 12))
        end
    else
        error("n must be 1 or 2")
    end
end
#alternative code for the Quark propagator
#=
#Self energy dressing function momentum integration Sigma_A
function SigmaA(A, M, p)
    global aold = zeros(Number, lastindex(p, 1))
    for k = 1:lastindex(p, 1)
        for j = 1:lastindex(A, 1)
            aold[k] +=  exp(-(sqrt(Complex(p[k])) - sqrt(u_q[j]))^2 * (eta / Gamma)^2) * (-2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * eta^2 * besselix(1, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / Gamma)^2) + (2*Gamma^2 + p[k] * eta^2 + u_q[j] * eta^2) * besselix(2, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / Gamma)^2)) / p[k] * 1 / (A[j] * (u_q[j] + M[j]^2)) *
                       (u_q[j]) * w_q[j]
        end
    end
    if lastindex(aold, 1) == 1
        return aold[1] .* eta^3/Gamma^2*Z_2^2
    end
    return aold .* eta^3/Gamma^2*Z_2^2
end
#Self energy dressing function momentum integration Sigma_M
function SigmaM(A, M, p)
    global aold = zeros(Number, lastindex(p, 1))
    for k = 1:lastindex(p, 1)
        for j = 1:lastindex(A, 1)
            aold[k] +=  exp(-(sqrt(Complex(p[k])) - sqrt(u_q[j]))^2 * (eta / Gamma)^2) * ((p[k] + u_q[j]) * besselix(1, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / Gamma)^2) - 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * besselix(2, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / Gamma)^2)) /sqrt(Complex(p[k])) * sqrt(u_q[j]) * M[j] / (A[j] * (u_q[j] + M[j]^2)) *
                        w_q[j]
        end
    end
    if lastindex(aold,1)==1
        return aold[1] .* eta^5/Gamma^2*Z_2^2
    end
    return aold .* eta^5/Gamma^2*Z_2^2
end
#Iteration scheme on the real axis grid until convergence withing epsilin<1E-9
function A_M(m,mu_s)
    global Z_2 = 0.97

    global A_q = ones(Float64, length(u_q))
    global A_old = ones(Number, length(u_q))
    global SigmaA_help = zeros(Number, length(u_q))
    global SigmaA_mu = 0

    global M_q = ones(Float64, length(u_q))
    global SigmaM_help = zeros(Number, length(u_q))
    global SigmaM_mu = 0
    global diff = 1

    while abs(diff) > 10^(-9)
        global A_old = A_q
        global SigmaA_help = SigmaA(A_q, M_q, u_q)
        global SigmaA_mu = SigmaA(A_q, M_q, mu_s)
        global Z_2 = 1 - SigmaA_mu
        global A_q = Z_2 .+ SigmaA_help 

        global diff = sqrt(sum(abs2.((A_old .- A_q))))

        global SigmaM_help = SigmaM(A_q, M_q, u_q)
        global SigmaM_mu = SigmaM(A_q, M_q, mu_s)
        global M_q = (m .+ SigmaM_help .- SigmaM_mu) ./ A_q
    end
    return A_q, M_q , Z_2, (m .- SigmaM_mu)
end

#Calculate the Quark Propagator dressing fucntion A & M
m_q=0.004 #current quark mass in GeV
mu_s=19 #renormalization potin in GeV^2
A1, M1, Z_2, RenormM = A_M(m_q, mu_s) #solution of the Dressing functions on the real grid

#Iterate 1 more time to calculate the dressing function on any complex momenta p

function A(p::Number)::Complex{Float64}
    return Z_2 .+  SigmaA(A1, M1, p) 
end
function M(p::Number)::Complex{Float64}
    return (RenormM .+  SigmaM(A1, M1, p)) ./ A2(p) 
end
=#

#Quark Propagator Fit provided by prof. Eichmann (for qualitative comparison)
function A2(p::Number)::Complex{Float64}
    return 0.95 + 0.3 / log((p) / (0.7)^2 + 2) + 0.1 / (1 + (p) / (0.7)^2) + 0.29 * exp(-0.1 * (p) / (0.7)^2) - 0.18 * exp(-3 * (p) / (0.7)^2)
end

function M2(p::Number)::Complex{Float64}
    return 0.06 / (1 + p / (0.7)^2) + 0.44 * exp(-0.66 * p / (0.7)^2) + 0.009 / (log(p / (0.7)^2 + 2))^(12 / 25)
end

#simplified effective interaction
function g(p::Number, q::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return Z_2^2 * 16 * pi / 3 * (pi * eta^7 / Gamma^4 * (sqrt(p)^2 + sqrt(q)^2 - 2*sqrt(p)*sqrt(q)*(z*zs + y*sqrt(1 - z^2)*sqrt(1 - zs^2))) * exp(-eta^2 * (sqrt(p)^2 + sqrt(q)^2 - 2*sqrt(p)*sqrt(q)*(z*zs + y*sqrt(1 - z^2)*sqrt(1 - zs^2))) / Gamma^2))
end
#----------------------------------------------------------------------------------
#include kernel entries from mathematica
#This needs to be adjusted depending on where your files are saved
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaK11.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaK22.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaK23.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaK32.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaK33.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaK44.jl")
#----------------------------------------------------------------------------------
#Everything for G Matrix
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG11.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG12.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG13.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG14.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG22.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG23.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG24.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG33.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG34.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliaG44.jl")
include("F:\\Julia\\2_Pion_BSE\\Expressions\\juliafactorG.jl")

#-----------------------------------------------------------------------------------------
#Calculating the Actual q^2 and z Integrals 
function a1(B_help1)
    aold = zeros(Number, lastindex(B_help1, 1), lastindex(B_help1, 2))
    for l = 1:lastindex(B_help1, 2)
        for k = 1:lastindex(B_help1, 1)
            for j = 1:lastindex(B_help1, 2)
                for i = 1:lastindex(B_help1, 1)
                    aold[k, l] += K_11[k,i,l,j]   * B_help1[i, j] * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold .* 0.5 ./ (2 .* pi) .^ 3
end
#for the updated eigenvalue at a specific point on the z-grid
function a1_lambda(B_help1,u_C_idx)
    aold = zeros(Number, lastindex(B_help1, 1))
    for k = 1:lastindex(B_help1, 1)
        for j = 1:lastindex(B_help1, 2)
            for i = 1:lastindex(B_help1, 1)
                aold[k] += K_11[k,i,u_C_idx,j] * B_help1[i, j] * (u_L1[i]) * w_L1[i] * w_C[j]
            end
        end
    end
    return aold .* 0.5 ./ (2 .* pi) .^ 3
end

function a2(B_help2, B_help3)
     aold = zeros(Number, lastindex(B_help2, 1), lastindex(B_help2, 2))
    for l = 1:lastindex(B_help2, 2)
        for k = 1:lastindex(B_help2, 1)
            for j = 1:lastindex(B_help2, 2)
                for i = 1:lastindex(B_help2, 1)
                    aold[k, l] += (K_22[k,i,l,j] * B_help2[i, j] + K_23[k,i,l,j] * B_help3[i, j]) * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold .* 0.5 ./ (2 .* pi) .^ 3
end

function a3(B_help2, B_help3)
    aold = zeros(Number, lastindex(B_help2, 1), lastindex(B_help2, 2))
    for l = 1:lastindex(B_help2, 2)
        for k = 1:lastindex(B_help2, 1)
            for j = 1:lastindex(B_help2, 2)
                for i = 1:lastindex(B_help2, 1)
                    aold[k, l] += (K_32[k,i,l,j] * B_help2[i, j] + K_33[k,i,l,j] * B_help3[i, j]) * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold .* 0.5 ./ (2 .* pi) .^ 3
end

function a4(B_help4)
    aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    aold[k, l] += K_44[k,i,l,j] * B_help4[i, j] * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold .* 0.5 ./ (2 .* pi) .^ 3
end
