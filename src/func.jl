# 
# This program defines the main functions  
# to be used in the program 
# Here we're using Nautral Units c = ħ = k_B = 1 
# Everything should be written in GeV
# 


# inporting the necessary libraries 
using QuadGK # this is the gaussian quadrature integral library


# defining the energy E_k 
function Ek(k,m=0.0)
    return sqrt(k^2+m^2)
end


# defininf the fermi momentum in terms of the 
# number density
function kfi(n::Float64,g::Float64=2.0)
    return cbrt(n*6*π^2/g)
    # return (6*π^2/g*n)^(1/3)
end

function ni(kf::Float64,g::Float64=2.0)
    return g/(6*π^2)*(kf^3)
end

# the chemical potential as a function of the fermi momentum
function μi(kf::Float64,m::Float64=0.0)
    return sqrt(kf^2 + m^2)
end


# defining the 3D integral of E_k
# or, the energy density of a free fermion gas
function εi(kf::Float64,M::Float64,g::Float64=2.0)
    return g*quadgk(k->k^2*Ek(k,M),0,kf)[1]/(2*π^2)
end


# defining the pressure 
function Piεi(n::Float64,m::Float64,g::Float64=2.0)
    kf = kfi(n,g)
    μ = μi(kf,m)
    ε = εi(kf,m,g)
    return [-ε+n*μ,ε]
end

function kfμ(μ::Float64,m::Float64)
    a = μ^2-m^2
    return real(sqrt(complex(a)))
end