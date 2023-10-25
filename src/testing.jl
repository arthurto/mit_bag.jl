# This program will be used to test 
# the routines that we are developing to
# solve the MIT Bag model.

begin
    using Plots 
    include("solve.jl")
end

begin
    Λ = 2
    initial = n_0*0.5

    x = LinRange(-Λ*1e-3,Λ*1e-3,100)
    y = LinRange(0,Λ,100)

    f = (x,y) -> eqs(initial,[y*n_0,x*n_0])[1]
    f1 = (x,y) -> eqs(initial,[y*n_0,x*n_0])[2]
    # f2 = (x,y) -> eqs(n_0,[y*n_0,x*n_0,1e-5*n_0])[3]
    z = @. f(x',y)
    z1 = @. f1(x',y)
    # z2 = @. f2(x',y)
    contour(x,y,z,levels = [0.0],color="black")
    contour!(x,y,z1,levels = [0.0],color="red")
    # contour!(x,y,z2,levels = [0.0],color="blue")
    # savefig("effective_mass.png")
end

begin
    x = LinRange(0,5,100)

    plot(x,x->thermodynamics(x*n_0)[2]/x)
    plot!(x,x->thermodynamics(x*n_0)[3]/x)
    plot!(x,x->thermodynamics(x*n_0)[4]/x)
    plot!(x,x->thermodynamics(x*n_0)[5]/x)
end

# 
# This program is used to test if the solutions are working properly
# 

begin
    
    Λ = 0.5
    N = 500
    x = [LinRange(5e-3,1e-1,N);LinRange(1e-1,Λ,N)]
    plot(x,(x->1/(3n_0)*thermodynamics(x*n_0)[2]/x).(x),
        label = "\$\\gamma_{u}\$",xlabel= "\$n_B/n_0\$",ylabel = "\$n_i/\\sum_i n_i\$",yaxis=:log,
        ylims=(1e-3,1e0),xlims=(0,Λ))
    plot!(x,(x->1/(3n_0)*thermodynamics(x*n_0)[3]/x).(x),
    label = "\$\\gamma_{d}\$")
    plot!(x,(x->1/(3n_0)*thermodynamics(x*n_0)[4]/x).(x),
    label = "\$\\gamma_{s}\$")
    plot!(x,(x->1/(3n_0)*100thermodynamics(x*n_0)[5]/x).(x),
    label = "\$100\\times \\gamma_{e}\$")
    # plot!(x,(x->1/(3n_0)*thermodynamics(x*n_0)[6]/x).(x),
    # label = "\$n_{\\mu}\$")
end
