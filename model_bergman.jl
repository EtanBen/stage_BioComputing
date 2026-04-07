#Bergman  


using Catalyst, Plots, OrdinaryDiffEq
@variables t 
function vd_dx(bw)
    """
        Compute the volume of distribution of D-Xylose (dL/kg) from the body weight (kg)
    """
    -0.0912608*bw + 10.9434
end

Insuline(x) = -1490.6 ./ ((-932.3 ./ (x .+ 3.8836)) .- x)

bergman_officiel = @reaction_network begin 

    a, G --> 0
    1,X+G --> X
    B0,0--> G
    k3, X --> 0
    Insuline(t)*P3, 0 --> X 
end 


CI = [:G => 30000/30, :X=> 0]

parameters = [:a=> 0.049,:B0=>4.42,:k3=>0.091,:P3=>0.0000896]   # a = -p1 b0 = p4 k3 = -p2

t_span=(0,360)

prob = ODEProblem(bergman_officiel,CI,t_span,parameters) 

sol_bergman_officiel = solve(prob)


gluco_plot=plot(sol_bergman_officiel,xlim=(0,360),ylim=(0,360),yticks=0:10:360,xlabel="time (min)",ylabel="Concentration du Glucose (mg/dL)")
