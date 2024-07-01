using Plots, SpecialFunctions

#The Expected Pion Mass should be around 0.1396GeV
#This needs to be adjusted depending on where u saved the file
include("f:/Julia/2_Pion_BSE/Pion_BSE_functions_latest.jl")
 
#parameters for the integration routing
num1 = 64
num2 = 16
num3 = 300
L = 10^3
u_q, w_q = Gauss_Legendre1(L^2, num3)
u_L1, w_L1 = Gauss_Legendre1(L^2, num1)
u_y, w_y = Gauss_Legendre(-1, 1, num2)
u_C, w_C = Gauss_Chebyshev2(-1, 1, num2)
#Paramters
global m = 0.004; 
global eta = 1.8 #1.92;
global Gamma = 0.72;
global mu_s= 19; #32
global Z_2=1 #initial value, 

#Routine for the Quark propagator
global A_help = map.(x -> 1, u_q)
global A_new= map.(x -> 1, u_q)
global M_help = map.(x -> 0.01, u_q)
global diff = 1 
while abs(diff) > 10^(-12) #breaking condition 
    global A_help=A_new
    global Sigma_A(x::Number)::Number =  eta^3/Gamma^2*Z_2^2 *
                    sum((exp.(-(sqrt.(Complex(x)) .- sqrt.(u_q)).^2 * (eta / Gamma)^2) .* (-2 * sqrt.(Complex(x)) .* sqrt.(u_q) * eta^2 .* besselix.(1, 2 * sqrt.(Complex(x)) .* sqrt.(u_q) * (eta / Gamma)^2) .+ (2*Gamma^2 .+ x * eta^2 .+ u_q * eta^2) .* besselix.(2, 2 * sqrt.(Complex(x)) .* sqrt.(u_q) * (eta / Gamma)^2)) ./ x .* 1 ./ (A_help .* (u_q .+ M_help.^2)) .*
                    u_q).*w_q)
   # global diff_A=sum(abs.(map.(A, u) .- A_help))               
    global Sigma_M(x::Number)::Number = eta^5/Gamma^2*Z_2^2*
                    sum((exp.(-(sqrt.(Complex(x)) .- sqrt.(u_q)).^2 * (eta / Gamma)^2) .* ((x .+ u_q) .* besselix.(1, 2 * sqrt.(Complex(x)) .* sqrt.(u_q) * (eta / Gamma)^2) .- 2 * sqrt.(Complex(x)) .* sqrt.(u_q) .* besselix.(2, 2 * sqrt.(Complex(x)) .* sqrt.(u_q) * (eta / Gamma)^2)) ./sqrt.(Complex(x)) .* M_help ./ (A_help .* (u_q .+ M_help.^2)) .*
                    sqrt.(u_q)) .* w_q)
  #  global diff_B=sum(abs.(map.(B_new, u) .- B_help))  
    global Z_2= 1 .- map(Sigma_A, mu_s)[1]
    global A_new = Z_2 .+ map.(Sigma_A, u_q)
    global M_help = (m .+  map.(Sigma_M, u_q) .- map.(Sigma_M, mu_s)[1]) ./ A_new
    global diff = sqrt(sum(abs2.((A_new .- A_help)))) # convergence check
end
A(x::Number)::Number = Z_2 + Sigma_A(x)
M(x::Number)::Number = (m .- map(Sigma_M, mu_s)[1] .+ Sigma_M(x))./A(x)
println("Quark Propagator is ready")

#perfoming the y integral first
function K_11s(p::Number, q::Number, z::Number, zs::Number)::Complex{Float64}
    return sum(g.(p, q, z, zs, u_y) .* K_11ss.(p, q, z, zs, u_y) .* w_y)
end

function K_22s(p::Number, q::Number, z::Number, zs::Number)::Complex{Float64}
    return sum(g.(p, q, z, zs, u_y) .* K_22ss.(p, q, z, zs, u_y) .* w_y)
end

function K_23s(p::Number, q::Number, z::Number, zs::Number)::Complex{Float64}
    return sum(g.(p, q, z, zs, u_y) .* K_23ss.(p, q, z, zs, u_y) .* w_y)
end

function K_32s(p::Number, q::Number, z::Number, zs::Number)::Complex{Float64}
    return sum(g.(p, q, z, zs, u_y) .* K_32ss.(p, q, z, zs, u_y) .* w_y)
end

function K_33s(p::Number, q::Number, z::Number, zs::Number)::Complex{Float64}
    return sum(g.(p, q, z, zs, u_y) .* K_33ss.(p, q, z, zs, u_y) .* w_y)
end

function K_44s(p::Number, q::Number, z::Number, zs::Number)::Complex{Float64}
    return sum(g.(p, q, z, zs, u_y) .* K_44ss.(p, q, z, zs, u_y) .* w_y)
end
#-----------------------------------------------------------------------------------------------------------------
#Creating our Grid for K Matrix for the values K(p,q,z,zs)=K(u_L1,u_L1,u_C,u_C)
#Creating a grid is faster than using the functions K_11s(p,q,z,zs)!
K_11 = Array{Number}(undef, lastindex(u_L1), lastindex(u_L1), lastindex(u_C), lastindex(u_C))
K_22 = Array{Number}(undef, lastindex(u_L1), lastindex(u_L1), lastindex(u_C), lastindex(u_C))
K_23 = Array{Number}(undef, lastindex(u_L1), lastindex(u_L1), lastindex(u_C), lastindex(u_C))
K_32 = Array{Number}(undef, lastindex(u_L1), lastindex(u_L1), lastindex(u_C), lastindex(u_C))
K_33 = Array{Number}(undef, lastindex(u_L1), lastindex(u_L1), lastindex(u_C), lastindex(u_C))
K_44 = Array{Number}(undef, lastindex(u_L1), lastindex(u_L1), lastindex(u_C), lastindex(u_C))
global idx_var = 0
for l = 1:lastindex(u_C)
    for k = 1:lastindex(u_L1)
        for j = 1:lastindex(u_C)
            for i = 1:lastindex(u_L1)
                K_11[k,i,l,j] = K_11s(u_L1[k], u_L1[i], u_C[l], u_C[j])    
                K_22[k,i,l,j] = K_22s(u_L1[k], u_L1[i], u_C[l], u_C[j])
                K_33[k,i,l,j] = K_33s(u_L1[k], u_L1[i], u_C[l], u_C[j])
                K_44[k,i,l,j] = K_44s(u_L1[k], u_L1[i], u_C[l], u_C[j])
                K_23[k,i,l,j] = K_23s(u_L1[k], u_L1[i], u_C[l], u_C[j])
                K_32[k,i,l,j] = K_32s(u_L1[k], u_L1[i], u_C[l], u_C[j])
            end
        end
    end
end
println("Kernel matrix is set up") #This takes up most of the time, there is alot of room left for efficiency of the code

#running the Integral for various values for the mass to see where lambda(m)=1
#here only the m=0.14 is written, but u could extend it
#note that only the last set of dressing functions are saved here, so you may modify the code
massspace = [0.14];
Eigenwert = []
for i in massspace
    global mass = i
   
    G_11 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_22 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_33 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_44 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_12 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_13 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_14 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_23 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_24 = Matrix{Number}(undef, length(u_L1), length(u_C))
    G_34 = Matrix{Number}(undef, length(u_L1), length(u_C))

    for j = 1:lastindex(u_C)
        for i = 1:lastindex(u_L1)
            G_11[i, j] = G_11s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_22[i, j] = G_22s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_33[i, j] = G_33s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_44[i, j] = G_44s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_12[i, j] = G_12s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_13[i, j] = G_13s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_14[i, j] = G_14s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_23[i, j] = G_23s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_24[i, j] = G_24s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
            G_34[i, j] = G_34s(u_L1[i], u_C[j])*factorG(u_L1[i], u_C[j])
        end
    end


    global A_help1 = ones(Number, length(u_L1), length(u_C))
    global A_help2 = ones(Number, length(u_L1), length(u_C))
    global A_help3 = ones(Number, length(u_L1), length(u_C))
    global A_help4 = ones(Number, length(u_L1), length(u_C))
    global lambda = 1

    for k = 1:10
        B_help1 = (G_11 .* A_help1 .+ G_12 .* A_help2 .+ G_13 .* A_help3 .+ G_14 .* A_help4)
        B_help2 = (G_12 .* A_help1 .+ G_22 .* A_help2 .+ G_23 .* A_help3 .+ G_24 .* A_help4)
        B_help3 = (G_13 .* A_help1 .+ G_23 .* A_help2 .+ G_33 .* A_help3 .+ G_34 .* A_help4)
        B_help4 = (G_14 .* A_help1 .+ G_24 .* A_help2 .+ G_34 .* A_help3 .+ G_44 .* A_help4)
        
        global lambda =  a1_lambda(B_help1,1)[1] / A_help1[1,1]
        println(lambda)
        
        global A_help1 = a1(B_help1) / lambda
        global A_help2 = a2(B_help2, B_help3) / lambda
        global A_help3 = a3(B_help2, B_help3) / lambda
        global A_help4 = a4(B_help4) / lambda

    end
    push!(Eigenwert, lambda)
end

#Calculate and plot the dressing functions f1,f2,f3,f4 in the original Tensor basis at the point z=u_C[9]=0.092
f1 = real.(A_help1[:, 9])
f2 = real.((A_help3[:, 9] * u_C[9] .- A_help2[:, 9] * sqrt(1 - u_C[9]^2)) / (mass * sqrt(1 - u_C[9]^2)))
f3 = -real.(A_help3[:, 9] ./ (mass * sqrt(1 - u_C[9]^2) * u_C[9] * u_L1))
f4 = real.(A_help4[:, 9] ./ (2 * im * mass * sqrt(1 - u_C[9]^2) * sqrt.(u_L1)))

p1 = scatter(u_L1, f1, xaxis=:log10, xlabel="\$p^2\$", label="\$f_1(p^2,z \\approx 0)\$")
scatter!(p1, u_L1, f2, xaxis=:log10, xlabel="\$p^2\$", label="\$f_2(p^2,z \\approx 0)\$")
scatter!(p1, u_L1, f3, xaxis=:log10, xlabel="\$p^2\$", label="\$f_3(p^2,z \\approx 0)\$")
scatter!(p1, u_L1, f4, xaxis=:log10, xlabel="\$p^2\$", label="\$f_4(p^2,z \\approx 0)\$")
xticks!([10^(-4), 10^(-3), 10^(-2), 10^(-1), 10^(0), 10^(1), 10^2, 10^3, 10^4, 10^5])
title!("Dressing function for m_{\\pi} = (" * string(mass) * ") GeV and eta = "* string(eta))


#p1 = plot(u_q, real.(A.(u_q)), xaxis=:log10, xlabel="\$p^2\$", label="\$M^2(p^2)\$")


#x_complex= collect((-10:0.01:0))
#p1 = plot(x_complex, abs2.(A.(x_complex)), xlabel="\$p^2\$", label="\$M^2(p^2)\$")
#xticks!(-10:1:0)

#ylims!(0,100)
#=
global f1_full = ones(Number, length(u_L1), length(u_C))
global f2_full = ones(Number, length(u_L1), length(u_C))
global f3_full = ones(Number, length(u_L1), length(u_C))
global f4_full = ones(Number, length(u_L1), length(u_C))
for i=1:lastindex(u_C)
    f1_full[:,i]= real.(A_help1[:, i])
    f2_full[:,i] = real.((A_help3[:, i] * u_C[i] .- A_help2[:, i] * sqrt(1 - u_C[i]^2)) / (mass * sqrt(1 - u_C[i]^2)))
    f3_full[:,i] = -real.(A_help3[:, i] ./ (mass * sqrt(1 - u_C[i]^2) * u_C[i] * u_L1))
    f4_full[:,i] = real.(A_help4[:, i] ./ (2 * im * mass * sqrt(1 - u_C[i]^2) * sqrt.(u_L1)))
end

#x_complex=collect(LinRange(-10,0,1000))
#p3=scatter(x_complex, abs2.(M(x_complex)), xlabel="\$p^2\$", label="\$M_{DSE}(p^2)\$" )
#scatter!(p3,x_complex,abs2.(M2.(Complex.(x_complex))), xlabel="\$p^2\$", label="\$M_{Fit}(p^2)\$")
#title!("Time-Like Quark Dressing Functions for \$\\eta = 1.8\$")

plotew = zeros(lastindex(Eigenwert))
for i = 1:lastindex(Eigenwert)
    plotew[i] = real(Eigenwert[i][1])
end

global Helpnorm = Vector{Number}(undef, length(u_L1) * length(u_C))
for j = 1:lastindex(u_C)
    for i = 1:lastindex(u_L1)
        idx_var1 = i + lastindex(u_L1) * (j - 1) 
        Helpnorm[idx_var1] = (16 *(4 *M(u_L1[i] - massspace[1]^2 / 4 - im * sqrt(u_L1[i]) * massspace[1] * u_C[j]) * M(u_L1[i] - massspace[1]^2 / 4 + im * sqrt(u_L1[i]) * massspace[1] * u_C[j]) *(f1_full[i,j]^2-f2_full[i,j]^2-f3_full[i,j]^2+f4_full[i,j]^2)+f1_full[i,j]^2 *(massspace[1]^2+4 * u_L1[i])-8 *im* f1_full[i,j]* f4_full[i,j]* massspace[1]* sqrt(u_L1[i]) * sqrt(1-u_C[j]^2)-f2_full[i,j]^2 *(massspace[1]^2+4 * u_L1[i] *(2* u_C[j]^2-1))-16 *f2_full[i,j]* f3_full[i,j]* u_L1[i]* u_C[j] *sqrt(1-u_C[j]^2)+f3_full[i,j]^2* massspace[1]^2+8* f3_full[i,j]^2 *u_L1[i]* u_C[j]^2-4* f3_full[i,j]^2 * u_L1[i]-f4_full[i,j]^2* massspace[1]^2-4 *f4_full[i,j]^2 * u_L1[i]))/(A(u_L1[i] - massspace[1]^2 / 4 - im * sqrt(u_L1[i]) * massspace[1] * u_C[j]) * A(u_L1[i] - massspace[1]^2 / 4 + im * sqrt(u_L1[i]) * massspace[1] * u_C[j])* (massspace[1]^2+4 *im* massspace[1]* sqrt(u_L1[i]) *u_C[j]-4* M(u_L1[i] - massspace[1]^2 / 4 - im * sqrt(u_L1[i]) * massspace[1] * u_C[j])^2-4* u_L1[i])* (massspace[1]^2-4* im* massspace[1]* sqrt(u_L1[i])* u_C[j]-4* M(u_L1[i] - massspace[1]^2 / 4 + im * sqrt(u_L1[i]) * massspace[1] * u_C[j])^2-4* u_L1[i]))
    end
end
global PionBSEnormalization= 0
for j = 1:lastindex(u_C)
    for i = 1:lastindex(u_L1)
        idx = i + lastindex(u_L1) * (j - 1) 
       global PionBSEnormalization += Helpnorm[idx]  *u_L1[i]* w_L1[i] * w_C[j]* 0.5 ./ (2 .* pi) .^ 3
    end
end
PionBSEnormalization *= 1.152510963349
BSEnormfactor=PionBSEnormalization^(-0.5)
global Helpdecay = Vector{Number}(undef, length(u_L1) * length(u_C))
for j = 1:lastindex(u_C)
    for i = 1:lastindex(u_L1)
        idx_var1 = i + lastindex(u_L1) * (j - 1) 
       global Helpdecay[idx_var1] = (16 * massspace[1]* (M(u_L1[i] - massspace[1]^2 / 4 - im * sqrt(u_L1[i]) * massspace[1] * u_C[j]) * (2 * im * f1_full[i,j] * massspace[1]+4 *f1_full[i,j]* sqrt(u_L1[i])* u_C[j]+4* im* f2_full[i,j]* M(u_L1[i] - massspace[1]^2 / 4 + im * sqrt(u_L1[i]) * massspace[1] * u_C[j])+4 * f4_full[i,j] *sqrt(u_L1[i])* sqrt(1-u_C[j]^2))+M(u_L1[i] - massspace[1]^2 / 4 + im * sqrt(u_L1[i]) * massspace[1] * u_C[j]) *(2 *im* f1_full[i,j]* massspace[1]-4 *f1_full[i,j]* sqrt(u_L1[i])* u_C[j]+4* f4_full[i,j]* sqrt(u_L1[i])* sqrt(1-u_C[j]^2))+im* (f2_full[i,j]* (massspace[1]^2+4* u_L1[i]* (2* u_C[j]^2-1))+8* f3_full[i,j]* u_L1[i]* u_C[j] *sqrt(1-u_C[j]^2))))/(A(u_L1[i] - massspace[1]^2 / 4 - im * sqrt(u_L1[i]) * massspace[1] * u_C[j]) *A(u_L1[i] - massspace[1]^2 / 4 + im * sqrt(u_L1[i]) * massspace[1] * u_C[j]) *(massspace[1]^2+4 *im* massspace[1]* sqrt(u_L1[i])* u_C[j]-4* M(u_L1[i] - massspace[1]^2 / 4 - im * sqrt(u_L1[i]) * massspace[1] * u_C[j])^2-4* u_L1[i])* (massspace[1]^2-4* im* massspace[1]* sqrt(u_L1[i])* u_C[j]-4* M(u_L1[i] - massspace[1]^2 / 4 + im * sqrt(u_L1[i]) * massspace[1] * u_C[j])^2-4* u_L1[i]))
    end
end
global PionBSEdecay= 0
for j = 1:lastindex(u_C)
    for i = 1:lastindex(u_L1)
        idx = i + lastindex(u_L1) * (j - 1) 
       global PionBSEdecay += Helpdecay[idx]  *u_L1[i]* w_L1[i] * w_C[j]* 0.5 ./ (2 .* pi) .^ 3
    end
end
PionBSEdecay *= im*Z_2*sqrt(8)/massspace[1]^2
#derivative of dlambda(m^2)/dm^2 = 1.152510963349 with the help of qti plot fitting the function
#lambda(m^2)=1,152510963349*m^2+0,9855125875343
#rawdata
#real.(-massspace .^ 2)
#9-element Vector{Float64}:
#-0.010000000000000002
#-0.0121
#-0.0144
#-0.016900000000000002
#-0.0225
#-0.0256
#-0.028900000000000006
#-0.0324
#-0.019600000000000003
#plotew
#9-element Vector{Float64}:
#0.9970946927346641
#0.9994805612535267
#1.002101798242948
#1.0049585835566677
#1.0114035209953403
#1.0149899567811995
#1.0188246682738782
#1.0229189975708801
#1.008058508109467
 p4 = scatter(real.(-massspace .^ 2), plotew, xlabel="\$P^2\$", ylabel="\$\\lambda\$")
scatter!(p4, [-0.1396^2], [1], label="Expected Value")
ylims!(0, 1.5)
xlims!(-1, 0.1)
#scatter!(p2, [-0.1396^2 / 4], [1])
#savefig("Eigenvalue1.pdf")
=#

#savefig("Dressing_functions_test2.pdf")
