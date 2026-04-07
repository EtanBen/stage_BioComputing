


using CSV,Optim,OrdinaryDiffEq,PEtab, DataFrames, Plots, Distributions, Statistics, Catalyst, Optimization,OptimizationNLopt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Import des datasets  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



dataXpostResec = CSV.read("./P_RT_X_postResecA.csv",DataFrame)
dataXpostResect    = dataXpostResec[2:end,1]
dataXpostResecmean = mean(Array(dataXpostResec[2:end,2:end]),dims = 2)
dataXpostResecmin  = minimum(Array(dataXpostResec[2:end,2:end]),dims = 2)
dataXpostResecmax  = maximum(Array(dataXpostResec[2:end,2:end]),dims = 2)
specimen = dataXpostResec[!,"GR16"][2:end]   # J'ajoute un seul specimen pour faire les premiers tests
BW = mean(skipmissing(Array(dataXpostResec[1,:])))





dataXPostJej =  CSV.read("./P_RT_X_postJejA.csv",DataFrame)
dataXPostJejTimes = dataXPostJej[2:end,1]
dataXPostJejMean = mean(Array(dataXPostJej[2:end,2:end]),dims = 2)








dataXPreResec = CSV.read("./P_RT_X_preResecA.csv",DataFrame)
dataXPreressecTimes   = dataXPreResec[2:end,1]
dataXPreResecMean = mean(Array(dataXPreResec[2:end,2:end]),dims = 2)

postRessecBWMean = mean(skipmissing(Array(dataXpostResec[1,:])))
preRessecBWMean = mean(skipmissing(Array(dataXPreResec[1,:])))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Model avec un seul compartiment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function vd_dx(bw)
    """
        Compute the volume of distribution of D-Xylose (dL/kg) from the body weight (kg)
    """
    -0.0912608*bw + 10.9434
end


model_1compartiment = @reaction_network begin 

     k_empt, Xs --> Xg
     alpha*k_abs,Xg --> Xp 
     k_elim, Xp --> 0 

end 

CI = [:Xs => 7500/(vd_dx(preRessecBWMean)*preRessecBWMean), :Xg => 0, :Xp => 0 ]
t_span = (0,360)
parameters =  [:k_empt => 0.0374, :k_abs => 0.222 , :k_elim => 0.00708, :alpha => 1] 

ode = ODEProblem(model_1compartiment, CI, t_span, parameters)

sol_1 = solve(ode) 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Model avec 5 compartiments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@variables t  


model_5compartiments = @reaction_network begin 

      k_empt, Xs --> Xg1      

      (a1)*k_abs, Xg1 --> Xp 
      k_trans, Xg1 --> Xg2 

      (a2)*k_abs, Xg2 --> Xp 
      k_trans, Xg2 --> Xg3

      (a3)*k_abs, Xg3 --> Xp 
      k_trans, Xg3 --> Xg4 

      (a4)*k_abs, Xg4 --> Xp 
      k_trans, Xg4 --> Xg5

      (a5)*k_abs, Xg5 --> Xp 
      

    
     

      k_elim, Xp --> 0
end 

CI_2 = [:Xs => 30000/(vd_dx(preRessecBWMean)*preRessecBWMean), :Xp => 0, :Xg1 =>0,:Xg2 => 0,:Xg3 => 0,:Xg4 => 0,:Xg5 => 0  ]  # Conditions initiales 

# coef "fonctionels" mais il y a  un alpha < 0 ( issue de l'optimiation de la julia doc )
alphaN = [ 0.28062839408874274, -0.11469799955980656, 0.062265161764182654, 0.40517098743019975, 0.31388702452191797]
coef_normalisation = sum((Array(alphaN)))
normalized_alphaN = (alphaN)./coef_normalisation



t_span = (0,500)

true_k = [0.0379,0.222,30/1100,0.00708] # k issue du papier 
parameters_k_estimed =[0.03067764735795039, 0.15175160666296494, 30/1100, 0.007081104488847296 ] # dans l'ordre k_empt  k_abs  k_trans  k_elim 

final_parameters = Pair.([:k_empt,:k_abs,:k_trans,:k_elim,:a1,:a2,:a3,:a4,:a5],[true_k;normalized_alphaN] )


# coef issures de l'optimisation avec PEtab.jl 
alpha_PE = [0.12200826304219005,0.001000000001743044,0.0010000000040480642, 0.02231382764493736,0.06016373826396602]
coef_NormPE = sum((Array(alpha_PE)))
alpha_NormPE = (alpha_PE)./coef_NormPE
parameters_PE  = Pair.([:k_empt,:k_abs,:k_trans,:k_elim,:a1,:a2,:a3,:a4,:a5],[true_k;alpha_PE] )  

prob_2 = ODEProblem(model_5compartiments,CI_2,t_span,parameters_PE)
sol_2 = solve(prob_2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   Determination des parametres   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# Premiere tentative d'optimiation avec le code de la doc Julia

function objective_function(p, _)
p = Pair.([:k_empt,:k_abs,:k_elim,:a1,:a2,:a3,:a4,:a5], p)
    opti_prob = remake(prob_2; p)
    sol = solve(opti_prob; saveat = times, save_idxs = :Xp, verbose = false, maxiters = 10000)
    SciMLBase.successful_retcode(sol) || return Inf
    return sum((sol .- dataXPreResecMean) .^2)
end

lbound = [0,0,0,0,0,0,0,0]
ubound = [1,1,1,1,1,1,1,1]

p_guess = [.1,.1,.1,.1,.1,.1,.1,.1]
pb_opti = OptimizationProblem(objective_function,p_guess,ub=ubound,lb=lbound)
sol_opti = solve(pb_opti,NLopt.LN_NELDERMEAD())   
"""
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   PEtab parameters fitting   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@unpack Xp = model_5compartiments 
obs_xp = PEtabObservable("obs_xp",Xp,0.5)
observables = Dict("obs_xp"=>obs_xp)
parameter_map=[:k_empt => 0.0374, :k_abs => 0.222 , :k_elim => 0.00708,:k_trans => 30/1100]  #parametre deja connu 

Pe_k_empt = PEtabParameter(:k_empt,lb = 0,ub =0.05)
Pe_k_abs = PEtabParameter(:k_abs,lb = 0,ub =0.5)
Pe_k_elim = PEtabParameter(:k_elim,lb = 0,ub =0.01)

Pe_a1 = PEtabParameter(:a1)
Pe_a2 = PEtabParameter(:a2 )
Pe_a3 = PEtabParameter(:a3 )
Pe_a4 = PEtabParameter(:a4 )
Pe_a5 = PEtabParameter(:a5)



params =[Pe_a1,Pe_a2,Pe_a3,Pe_a4,Pe_a5]
measurements = DataFrame(obs_id = "obs_xp", time = vec(dataXPreressecTimes), measurement = vec(dataXPreResecMean))
petab_model = PEtabModel(model_5compartiments, obs_xp, measurements, params  ,speciemap  = CI_2,parametermap = parameter_map)
petab_problem = PEtabODEProblem(petab_model)
p0 = [

    0.1,    
    0.1,   
    0.1,    
    0.1,    
    0.1     
]

res = calibrate(petab_problem, p0, IPNewton())



#print(get_ps(res, petab_problem))   #Pour imprimer les alpha optimise 






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Graph des resulats     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




plot(sol_2.t,sol_2[:Xp], label="\$model 5 compartiments~X_p\$", ylabel = "\$ X_p (mg/dL)\$", xlabel ="\$times~(min)\$")
#plot!(times,specimen,label="specimen",m=:star,l=:dash)
plot!(dataXPreressecTimes,dataXPreResecMean,m=:circle,l=:dash,label="\$Data~avant~opération\$")
#plot!(sol_1.t,sol_1[:Xp],label="Model 1 compartiment")
#plot!(times,dataXpostResecmean,ylims = (0.0,55.0),c = 4,m = :circle,l=:dash,ribbon = (dataXpostResecmean-dataXpostResecmin,dataXpostResecmax-dataXpostResecmean),label = "\$data~X_p\$")
#println("les parametres k opti sont ", sol_opti.u)


plot(sol_2.t,sol_2[:Xp])

#plot!(dataXPostJejTimes,dataXPostJejMean,m = :square)

