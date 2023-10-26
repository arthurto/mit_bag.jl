# 
# This program will generate the EOS as save it as a csv file 
# to be used to solve the TOV equations
#  

#
# Including the solving functions 
# and the necessary libraries
# 

begin
    include("solve.jl")
    using DataFrames
    using CSV 
end

# Defining the value of the bag constant to be 
# the one that the E/A ≈ 928 MeV
begin
    B = 0.1591^4 # GeV
end

# generating the data and saving it
begin
    function EosGenerate(;N=100,B=0.1591^4,max_dens = 15.0,min_dens=2.0)
        matrix = zeros(11,N) # matrix that stores the data to be written in the 
        nnb = LinRange(min_dens*n_0,max_dens*n_0,N)
        for i in 1:N 

            nb = nnb[i]

            matrix[:,i] = thermodynamics(nb,B)
        end 
        return matrix 
    end  
end

begin
    function main()
        matrix = EosGenerate()
        data = DataFrame(
            nb=matrix[1,:],
            nu=matrix[2,:],
            nd=matrix[3,:],
            ns=matrix[4,:],
            ne=matrix[5,:],
            nm=matrix[6,:],
            μu=matrix[7,:],
            μd=matrix[8,:],
            μe=matrix[9,:],
            P=matrix[10,:],
            e=matrix[11,:])
    
        CSV.write("eos_files/eos.csv",data)        
    end

    if abspath(PROGRAM_FILE) == @__FILE__
        main()
    end    
end
