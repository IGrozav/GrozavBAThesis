#global eta = 1.8
global p_0 = 0.72
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

#Quark Propagator Dressing Function selfmade
function bessel_map(n::Number, x::Number)::Number
    if real(x) < 0
        x = -x
    end
    if (n == 1)
        if real(x) >= 3.75
            return Complex(x) .^ (-0.5) .* (0.39894228 .- 0.03988024 .* (x ./ 3.75) .^ (-1) .- 0.00362018 .* (x ./ 3.75) .^ (-2) .+ 0.00163801 .* (x ./ 3.75) .^ (-3) .- 0.01031555 .* (x ./ 3.75) .^ (-4) .+ 0.02282967 .* (x ./ 3.75) .^ (-5) .- 0.02895312 .* (x ./ 3.75) .^ (-6) .+ 0.01787654 .* (x ./ 3.75) .^ (-7) .- 0.00420059 .* (x ./ 3.75) .^ (-8))
        else
            return exp.(-real.(x)) .* x .* (0.5 .+ 0.87890594 .* (x ./ 3.75) .^ 2 .+ 0.51498869 .* (x ./ 3.75) .^ 4 .+ 0.15084934 .* (x ./ 3.75) .^ 6 .+ 0.02658733 .* (x ./ 3.75) .^ 8 .+ 0.00301532 .* (x ./ 3.75) .^ 10 .+ 0.00032411 .* (x ./ 3.75) .^ 12)
        end
    end
    if (n == 2)
        if real(x) >= 3.75
            return Complex(x) .^ (-0.5) .* (0.39894228 .+ 0.01328582 .* (x ./ 3.75) .^ (-1) .+ 0.00225319 .* (x ./ 3.75) .^ (-2) .- 0.00157565 .* (x ./ 3.75) .^ (-3) .+ 0.00916281 .* (x ./ 3.75) .^ (-4) .- 0.02057706 .* (x ./ 3.75) .^ (-5) .+ 0.02635537 .* (x ./ 3.75) .^ (-6) .- 0.01647633 .* (x ./ 3.75) .^ (-7) .+ 0.00392377 .* (x ./ 3.75) .^ (-8)) .- 2 .* Complex(x) .^ (-1.5) .* (0.39894228 .- 0.03988024 .* (x ./ 3.75) .^ (-1) .- 0.00362018 .* (x ./ 3.75) .^ (-2) .+ 0.00163801 .* (x ./ 3.75) .^ (-3) .- 0.01031555 .* (x ./ 3.75) .^ (-4) .+ 0.02282967 .* (x ./ 3.75) .^ (-5) .- 0.02895312 .* (x ./ 3.75) .^ (-6) .+ 0.01787654 .* (x ./ 3.75) .^ (-7) .- 0.00420059 .* (x ./ 3.75) .^ (-8))
        else
            return exp.(-real.(x)) .* ((1 .+ 3.5156229 .* (x ./ 3.75) .^ 2 .+ 3.0899424 .* (x ./ 3.75) .^ 4 .+ 1.2067492 .* (x ./ 3.75) .^ 6 .+ 0.2659732 .* (x ./ 3.75) .^ 8 .+ 0.0360768 .* (x ./ 3.75) .^ 10 .+ 0.0045813 .* (x ./ 3.75) .^ 12) .- 2 .* (1 ./ 2 .+ 0.87890594 .* (x ./ 3.75) .^ 2 .+ 0.51498869 .* (x ./ 3.75) .^ 4 .+ 0.15084934 .* (x ./ 3.75) .^ 6 .+ 0.02658733 .* (x ./ 3.75) .^ 8 .+ 0.00301532 .* (x ./ 3.75) .^ 10 .+ 0.00032411 .* (x ./ 3.75) .^ 12))
        end
    else
        error("n must be 1 or 2")
    end
end
#=
function SigmaA(A, M, p)
    global aold = zeros(Number, lastindex(p, 1))
    for k = 1:lastindex(p, 1)
        for j = 1:lastindex(A, 1)
            aold[k] += 3 * 0.72^2* pi / (2* eta^4) * exp(-(sqrt(Complex(p[k])) - sqrt(u_q[j]))^2 * (eta / 0.72)^2) * (-2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * eta^2 * bessel_map(1, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / 0.72)^2) + (2*0.72^2 + p[k] * eta^2 + u_q[j] * eta^2) * bessel_map(2, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / 0.72)^2)) / p[k] * 1 / (A[j] * (u_q[j] + M[j]^2)) *
                       (u_q[j]) * w_q[j]
        end
    end
    if lastindex(aold, 1) == 1
        return aold[1] .* Z_2^2 .* 16 .* pi^2 .* eta .^ 7 ./ (3 .* 0.72 .^ 4 .* (2 .* pi) .^ 3)
    end
    return aold .* Z_2^2 .* 16 .* pi^2 .* eta .^ 7 ./ (3 .* 0.72 .^ 4 .* (2 .* pi) .^ 3)
end

function SigmaM(A, M, p)
    global aold = zeros(Number, lastindex(p, 1))
    for k = 1:lastindex(p, 1)
        for j = 1:lastindex(A, 1)
            aold[k] += pi / 2 * (0.72 / eta)^2 * exp(-(sqrt(Complex(p[k])) - sqrt(u_q[j]))^2 * (eta / 0.72)^2) * ((p[k] + u_q[j]) * bessel_map(1, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / 0.72)^2) - 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * bessel_map(2, 2 * sqrt(Complex(p[k])) * sqrt(u_q[j]) * (eta / 0.72)^2)) / ( sqrt(Complex(p[k])) * sqrt(u_q[j]) ) * M[j] / (A[j] * (u_q[j] + M[j]^2)) *
                       (u_q[j]) * w_q[j]
        end
    end
    if lastindex(aold,1)==1
        return aold[1] .* 3 .* Z_2^2 .* 16 .* pi^2 .* eta .^ 7 ./ (3 .* 0.72 .^ 4 .* (2 .* pi) .^ 3)
    end
    return aold .* 3 .* Z_2^2 .* 16 .* pi^2 .* eta .^ 7 ./ (3 .* 0.72 .^ 4 .* (2 .* pi) .^ 3)
end

function A_M(m, mu_s)
    global Z_2 = 0.97

    global A_q = ones(Float64, length(u_q))
    global A_old = ones(Number, length(u_q))
    global SigmaA_help = zeros(Number, length(u_q))
    global SigmaA_mu = 0

    global M_q = ones(Float64, length(u_q))
    global SigmaM_help = zeros(Number, length(u_q))
    global SigmaM_mu = 0
    global diff = 1

    while real(diff) > 10^(-8)
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
#Quark Propagator Eichmann Fit
function A2(p::Number)::Complex{Float64}
    return 0.95 + 0.3 / log((p) / (0.7)^2 + 2) + 0.1 / (1 + (p) / (0.7)^2) + 0.29 * exp(-0.1 * (p) / (0.7)^2) - 0.18 * exp(-3 * (p) / (0.7)^2)
end

function M2(p::Number)::Complex{Float64}
    return 0.06 / (1 + p / (0.7)^2) + 0.44 * exp(-0.66 * p / (0.7)^2) + 0.009 / (log(p / (0.7)^2 + 2))^(12 / 25)
end

#Everything for K Matrix
function ls(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return x + xs - 2 * sqrt(x) * sqrt(xs) * (z * zs + y * sqrt(1 - z^2) * sqrt(1 - zs^2))
end

#k^2=p^2+q^2-2*p*q*(z*zs+y*sqrt(1-z^2)*sqrt(1-zs^2))
function g(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return Z_2^2 * 16 * pi / 3 * (pi * eta^7 * ls(x, xs, z, zs, y) / p_0^4 * exp(-eta^2 * ls(x, xs, z, zs, y) / p_0^2)) #+
                                 # 2 * pi * 12 / 25 * (1 - exp(-ls(x, xs, z, zs, y))) / ( ls(x, xs, z, zs, y) * log(exp(2) - 1 + (1 + ls(x, xs, z, zs, y) / 0.234^2)^2))
end

function K_11ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return -((-2 *sqrt(x)* sqrt(xs)* (2 *y^3* sqrt(1-z^2) *sqrt(1-zs^2)+(y^2+1) *z *zs)+x* (y^2+1)+xs* y^2+xs)/(2 *ls(x, xs, z, zs, y)))
end

function K_22ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (-2 *sqrt(x) *sqrt(xs) *((y^2+1)* z *zs-2 *y* sqrt(1-z^2) *sqrt(1-zs^2))+x *(y^2+1) *(2* z^2-1)+xs* (y^2+1)* (2 *zs^2-1))/(2 *ls(x, xs, z, zs, y))
end

function K_33ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (-sqrt(x)* sqrt(xs)* (y^2 *sqrt(1-z^2) *sqrt(1-zs^2)-2* y* z* zs+sqrt(1-z^2)* sqrt(1-zs^2))+x* (y-2 *y *z^2)+xs* y *(1-2 *zs^2))/ls(x, xs, z, zs, y)
end

function K_44ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (sqrt(x)* sqrt(xs) *(-3 *y^2 *sqrt(1-z^2)* sqrt(1-zs^2)-2* y* z* zs+sqrt(1-z^2) *sqrt(1-zs^2))+x* y+xs* y)/ls(x, xs, z, zs, y)
end

function K_55ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return  3*y
end

function K_66ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return y *(-((2 *(sqrt(x)* z-sqrt(xs) *zs)^2)/ls(x, xs, z, zs, y))-1)
end

function K_77ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (y* (2* sqrt(x)* sqrt(xs) *(2* y^2* sqrt(1-z^2)* sqrt(1-zs^2)+y* z *zs+sqrt(1-z^2) *sqrt(1-zs^2))+x* y *(2 *z^2-3)+xs* y *(2 *zs^2-3)))/ls(x, xs, z, zs, y)
end

function K_88ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (y *(-2 *sqrt(x)* sqrt(xs)* (y *z *zs+sqrt(1-z^2) *sqrt(1-zs^2))+x* y+xs *y))/ls(x, xs, z, zs, y)
end

function K_16ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (sqrt(xs)* (y^2-1) *sqrt(2-2*zs^2) *(sqrt(xs)* zs-sqrt(x) *z))/ls(x, xs, z, zs, y)
end

function K_61ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (sqrt(x)* (y^2-1) *sqrt(2-2*z^2)* (sqrt(x) *z-sqrt(xs) *zs))/ls(x, xs, z, zs, y)
end

function K_17ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return ((y^2-1)* (-2 *sqrt(x)* sqrt(xs)* (2 *y *sqrt(1-z^2) *sqrt(1-zs^2)+z* zs)+x+xs *(3-2* zs^2)))/(sqrt(2) * ls(x, xs, z, zs, y))
end

function K_71ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return -(((y^2-1) *(2 *sqrt(x)* sqrt(xs)* (2* y *sqrt(1-z^2)* sqrt(1-zs^2)+z* zs)+x* (2* z^2-3)-xs))/(sqrt(2) * ls(x, xs, z, zs, y)))
end

function K_23ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return ((sqrt(x)* z-sqrt(xs)* zs) *(2 *sqrt(x)* y *sqrt(1-z^2)-sqrt(xs)* (y^2+1) *sqrt(1-zs^2)))/ls(x, xs, z, zs, y)
end

function K_32ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return ((sqrt(x)* z-sqrt(xs) *zs) *(sqrt(x) *(y^2+1) *sqrt(1-z^2)-2 *sqrt(xs)* y* sqrt(1-zs^2)))/ls(x, xs, z, zs, y)
end

function K_76ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return (2* y *(sqrt(x)* z-sqrt(xs)* zs) *(sqrt(xs)* y* sqrt(1-zs^2)-sqrt(x)* sqrt(1-z^2)))/ls(x, xs, z, zs, y)
end

function K_67ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return -((2 *y *(sqrt(x) *z-sqrt(xs) *zs) *(sqrt(x)* y *sqrt(1-z^2)-sqrt(xs) *sqrt(1-zs^2)))/ls(x, xs, z, zs, y))
end

function K_28ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return K_71ss(x,xs,z,zs,y) + sqrt(2)*(1-y^2)
end

function K_82ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return K_17ss(x,xs,z,zs,y) + sqrt(2)*(1-y^2)
end

function K_38ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return -K_61ss(x,xs,z,zs,y)
end

function K_83ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return -K_16ss(x,xs,z,zs,y)
end

function K_99ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return 3
end

function K_1010ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return K_66ss(x,xs,z,zs,y)/y
end

function K_1011ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return K_67ss(x,xs,z,zs,y)/y
end

function K_1110ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return K_76ss(x,xs,z,zs,y)/y
end

function K_1111ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return K_77ss(x,xs,z,zs,y)/y
end

function K_1212ss(x::Number, xs::Number, z::Number, zs::Number, y::Number)::Complex{Float64}
    return K_88ss(x,xs,z,zs,y)/y
end
#----------------------------------------------------------------------------------
#Everything for G Matrix

function factor_G(x::Number, z::Number, Q::Number)::Complex{Float64}
    return A(x + Q^2/4 - sqrt(x)*Q*z) *(x + Q^2/4 - sqrt(x)*Q*z + (M(x + Q^2/4 - sqrt(x)*Q*z))^2) * A(x + Q^2/4 + sqrt(x)*Q*z) * (x + Q^2/4 + sqrt(x)*Q*z + (M(x + Q^2/4 + sqrt(x)*Q*z))^2) 
end

function G_11s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (x+M(x + Q^2/4 - sqrt(x)*Q*z) * M(x + Q^2/4 + sqrt(x)*Q*z)-Q^2/4) / factor_G(x, z, Q)
end

function G_22s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (x* (2 *z^2-1)+M(x + Q^2/4 - sqrt(x)*Q*z) * M(x + Q^2/4 + sqrt(x)*Q*z)-Q^2/4) / factor_G(x, z, Q)
end

function G_33s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (x* (1-2 *z^2)+M(x + Q^2/4 - sqrt(x)*Q*z)* M(x + Q^2/4 + sqrt(x)*Q*z)+Q^2/4) / factor_G(x, z, Q)
end

function G_44s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (-x+M(x + Q^2/4 - sqrt(x)*Q*z) * M(x + Q^2/4 + sqrt(x)*Q*z)+Q^2/4) / factor_G(x, z, Q)
end

function G_12s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (1/2 * im * (M(x + Q^2/4 - sqrt(x)*Q*z) *(2 * sqrt(x)* z+Q) + M(x + Q^2/4 + sqrt(x)*Q*z) * (Q-2 * sqrt(x)* z))) / factor_G(x, z, Q)
end

function G_13s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return  (im * sqrt(x)* sqrt(1-z^2)* (M(x + Q^2/4 - sqrt(x)*Q*z)-M(x + Q^2/4 + sqrt(x)*Q*z))) / factor_G(x, z, Q)
end

function G_14s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (-sqrt(x)* Q *sqrt(1-z^2)) / factor_G(x, z, Q)
end

function G_23s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (2 * x* z * sqrt(1-z^2)) / factor_G(x, z, Q)
end

function G_24s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (im* sqrt(x)* sqrt(1-z^2)* (M(x + Q^2/4 - sqrt(x)*Q*z) + M(x + Q^2/4 + sqrt(x)*Q*z))) / factor_G(x, z, Q)
end

function G_34s(x::Number, z::Number, Q::Number)::Complex{Float64}
    return (-1/2 * im* (M(x + Q^2/4 - sqrt(x)*Q*z) * (2* sqrt(x)* z+Q) - M(x + Q^2/4 + sqrt(x)*Q*z) *(Q-2 *sqrt(x)* z))) / factor_G(x, z, Q)
end

#-----------------------------------------------------------------------------------------
#Calculating the Actual Integral
function a1(B_help1, B_help6, B_help7)
    global aold = zeros(Number, lastindex(B_help1, 1), lastindex(B_help1, 2))
    for l = 1:lastindex(B_help1, 2)
        for k = 1:lastindex(B_help1, 1)
            for j = 1:lastindex(B_help1, 2)
                for i = 1:lastindex(B_help1, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += (K_11[idx] * B_help1[i, j]  + K_16[idx] * B_help6[i, j] + K_17[idx] * B_help7[i, j])* (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a2(B_help2, B_help3, B_help8)
    global aold = zeros(Number, lastindex(B_help2, 1), lastindex(B_help2, 2))
    for l = 1:lastindex(B_help2, 2)
        for k = 1:lastindex(B_help2, 1)
            for j = 1:lastindex(B_help2, 2)
                for i = 1:lastindex(B_help2, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += (K_22[idx] * B_help2[i, j] + K_23[idx] * B_help3[i, j] + K_28[idx] * B_help8[i, j]) * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a3(B_help2, B_help3, B_help8)
    global aold = zeros(Number, lastindex(B_help2, 1), lastindex(B_help2, 2))
    for l = 1:lastindex(B_help2, 2)
        for k = 1:lastindex(B_help2, 1)
            for j = 1:lastindex(B_help2, 2)
                for i = 1:lastindex(B_help2, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += (K_32[idx] * B_help2[i, j] + K_33[idx] * B_help3[i, j] + K_38[idx] * B_help8[i, j]) * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a4(B_help4)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help4, 1) * (j - 1) + lastindex(B_help4, 1) * lastindex(B_help4, 2) * (k - 1) + lastindex(B_help4, 1)^2 * lastindex(B_help4, 2) * (l - 1)
                    aold[k, l] += (K_44[idx] * B_help4[i, j])* (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a5(B_help5)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += (K_55[idx] * B_help5[i, j])* (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a6(B_help1, B_help6, B_help7)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += ( K_61[idx] * B_help1[i, j]  +  K_66[idx] * B_help6[i, j] + K_67[idx] * B_help7[i, j]  )* (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a7(B_help1, B_help6, B_help7)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += ( K_71[idx] * B_help1[i, j]  + K_76[idx] * B_help6[i, j] + K_77[idx] * B_help7[i, j] )  * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a8(B_help2, B_help3, B_help8)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += (K_82[idx] * B_help2[i, j] + K_83[idx] * B_help3[i, j] + K_88[idx] * B_help8[i, j] )* (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a9(B_help9)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += K_99[idx] * B_help9[i, j] * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a10(B_help10, B_help11)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += (K_1010[idx] * B_help10[i, j] + K_1011[idx] * B_help11[i, j])* (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a11(B_help10,B_help11)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += (K_1110[idx] * B_help10[i, j] + K_1111[idx] * B_help11[i, j])* (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end

function a12(B_help12)
    global aold = zeros(Number, lastindex(B_help4, 1), lastindex(B_help4, 2))
    for l = 1:lastindex(B_help4, 2)
        for k = 1:lastindex(B_help4, 1)
            for j = 1:lastindex(B_help4, 2)
                for i = 1:lastindex(B_help4, 1)
                    idx = i + lastindex(B_help1, 1) * (j - 1) + lastindex(B_help1, 1) * lastindex(B_help1, 2) * (k - 1) + lastindex(B_help1, 1)^2 * lastindex(B_help1, 2) * (l - 1)
                    aold[k, l] += K_1212[idx] * B_help12[i, j] * (u_L1[i]) * w_L1[i] * w_C[j]
                end
            end
        end
    end
    return aold 
end
