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

# generating the data and saving it
begin
    N = 100 # number of points
    matrix = zeros(11,N) # matrix that stores the data to be written in the 
    # data frame 
    for i in 1:N 
        nb = 5n_0*i/N
        matrix[:,i] = thermodynamics(nb)
    end 
end

begin
    # Defining the dataframe that 
    # will be used to save the data 
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
