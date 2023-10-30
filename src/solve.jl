#
# This program will solve the equations for the 
# charge neutrality and the chemical equilibrium 
# of the matter composed of free fermions
# This program uses Nautral Units c = ħ = k_B = 1 
# Everything should be written in GeV
# 

# Importing the nonlinear solving package 
using NLsolve 

# Including 'func.jl'
include("func.jl")

# Constants 
# n_0 
n_0 = 0.15*(197^3*1e-9) # Nuclear saturation density in GeV^3

# Function that returns μu,μe,μd,nd,ns,nm
# as functions of nu and ne
function qts(nu,ne)

    # finding the chemical potential as a function of 
    # n for q_u and e^-
    μu = μi(kfi(nu,6.0),0.002) # m_u ≈ 0.002 [GeV]   g = (spin)×(colors)=6
    μe = kfi(ne,2.0) # m_μ ≈ 0.000 [GeV]   g = 2

    # using chemical equillibrium to
    # define μ_d and μ_s 
    μd = (μu + μe)

    nm = ni(kfμ(μe,0.105))     # since μ_μ = μ_e
    ns = ni(kfμ(μd,0.093),6.0) # since μs = μd 
    nd = ni(kfμ(μd,0.005),6.0) 


    return μu,μe,μd,nd,ns,nm
end

# defining the vector function to be solved
# vec(F(vec(x))) = vec(0)
# in this case we have only two equations to be solved
function eqs(nb,x) #
    
    # Getting the values for n_u and n_e 
    nu,ne = x
    
    # using only the necessary quantities
    nd,ns,nm = qts(nu,ne)[4:end]

    return [
            nb - (nu+nd+ns)/3,            # equation for defining the baryon density
            -nm - ne + (2*nu -nd -ns)/3,  # equation for charge neutrality
            ]
end

# Defining the solver function as a function of nb
function solv(nb)
    return nlsolve(x->eqs(nb,x),[0.5,1e-2]*nb).zero
end

# Defining the function that returns all the 
# thermodynamical quantities as a function of nb
function thermodynamics(nb,B = (0.154)^4)
    nu,ne = solv(nb)

    μu,μe,μd,nd,ns,nm = qts(nu,ne)

    Pu,εu = Piεi(nu,0.002,6.0) 
    Pd,εd = Piεi(nd,0.005,6.0) 
    Ps,εs = Piεi(ns,0.093,6.0) 
    Pe,εe = Piεi(ne,0.000) 
    Pm,εm = Piεi(nm,0.105)
    return [# 1, 2, 3, 4, 5, 6,
             nb,nu,nd,ns,ne,nm,
             #7, 8, 9,
             μu,μd,μe,
             # 10 = P 
             sum([Pu,Pd,Ps,Pe,Pm])-B, # Adding and subtracting
             # 11 = ε
             sum([εu,εd,εs,εe,εm])+B  # the bag constant 
            ] 
end