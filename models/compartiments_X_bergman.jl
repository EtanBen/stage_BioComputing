using CSV,Optim,OrdinaryDiffEq,PEtab, DataFrames, Plots, Distributions, Statistics, Catalyst, Optimization,OptimizationNLopt
include("visualiation_data.jl")
include("../dalla_man_CRN.jl")
using .itp_insuline
#@register_symbolic Insuline(t)
@variables t 
BW = 78 
insuline_dallaMan(t) = sol_dallaman(t; idxs=I_p)     #Insuline fournit par le model de DALLAMAN  
glucose_dallaMan(t) = sol_dallaman(t;idxs = G_p ) # Glucose fournit par le model de DALLAMAN



@register_symbolic insuline_dallaMan(t)
@register_symbolic glucose_dallaMan(t)


# volume de distribution 


function vd_g(bw)
    """
        Compute the volume of distribution of Glucose (dL/kg) from the body weight (kg)
        linear
    """
    -0.00899437*bw + 1.64445
end


function vd_dx(bw)
    """
        Compute the volume of distribution of D-Xylose (dL/kg) from the body weight (kg)
    """
    -0.0912608*bw + 10.9434
end
#Insuline2(x) = -1490.6 ./ ((-932.3 ./ (x .+ 3.8836)) .- x)


# glucose des cochons
GlucoData = CSV.read("./P_RT_G_preResecA.csv",DataFrame)
GlucoDataTime = GlucoData[2:end,1]
GlucoDataMean = mean(Array(GlucoData[2:end,2:end]),dims = 2)
BW_pig = mean(skipmissing(Array(GlucoData[1,:])))




val_gluco_dallaMan = glucose_dallaMan.(GlucoDataTime)./1.88  # on convertit les mg/kg en mg/dl via le Vg de dallaman !!! ATTENTION !!! il change si c'est un diab 

#println(val_gluco_dallaMan)




model_mix5comp = @reaction_network begin   #TODO changer les x par des g 

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
      

      1,X+Xp --> X
      ins_elim,X --> 0
      insuline_dallaMan(t)*ins_app, 0--> X
      k_elim, Xp --> 0
end 


# Bric a Brac 
mixed_parameters_PE =  [:a2 => 13.624954725940261, :ins_elim => 999.9999998934837, :a5 => 7.004193422048467, :a4 => 0.0010000000009834436, :k_trans => 0.02727272727272727, :ins_app => 0.0010000000001065203, :k_empt => 1.0329746820572028, :k_abs => 0.0010000000009834436, :B0 => 15.174879848678678, :a3 => 0.4502304160042123, :a1 => 2.3493343173824353, :k_elim => 0.19990528931742133]
coef_NormPE = sum((Array(alpha_PE)))
alpha_NormPE = (alpha_PE)./coef_NormPE
alpaha_02 = 0.2.*ones(5,1)




CI = [:Xs => 90000/(vd_g(BW)*BW), :Xp => 90, :X =>0,:Xg1 =>0,:Xg2 => 0,:Xg3 => 0,:Xg4 => 0,:Xg5 => 0  ]  # Conditions initiales   #  le xp init ??? 
true_k = [0.222,0.0379,30/1100,0.00628] # k issue du papier 
alpha_PE = [0.02215777201700488,0.002556170177838383,0.009999414260685416,0.0010000005301935998,0.0024545573258758203]   # alpha optimisation PETab + l'ancien k_abs ( 0.222)
parameters_old  = Pair.([:k_abs,:k_empt,:k_trans,:k_elim,:a1,:a2,:a3,:a4,:a5,:ins_app,:ins_elim ],[true_k;alpha_PE;[0.0000896,0.091]] )  # params insuline papier de bergmann  
PE_params = [:a2 => 0.0010000414002961195, :a5 => 0.0037041084449406366, :a4 => 0.0010000308476614068, :ins_elim => 0.091, :k_trans => 0.02727272727272727, :ins_app => 8.96e-5, :k_empt => 0.0379, :k_abs => 0.222, :a3 => 0.016702739932300845, :a1 => 0.03175557113587718, :k_elim => 0.00628]


t_span=(0,360)
prob_mix = ODEProblem(model_mix5comp,CI,t_span,PE_params)
sol_mix = solve(prob_mix)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  recherche de parametres via PEtab  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
@unpack Xp = model_mix5comp 
obs_xp = PEtabObservable("obs_xp",Xp,0.5)
observables = Dict("obs_xp"=>obs_xp)

parameter_map = [ :k_trans => 30/1100, :k_empt => 0.0379 , :k_elim => 0.00628, :ins_app => 0.0000896, :ins_elim => 0.091, :k_abs=>0.222]


Pe_k_empt = PEtabParameter(:k_empt)
Pe_k_abs = PEtabParameter(:k_abs)
Pe_k_elim = PEtabParameter(:k_elim)

Pe_a1 = PEtabParameter(:a1)
Pe_a2 = PEtabParameter(:a2 )
Pe_a3 = PEtabParameter(:a3 )
Pe_a4 = PEtabParameter(:a4 )
Pe_a5 = PEtabParameter(:a5)

#Pe_B0 = PEtabParameter(:B0)
Pe_ins_app = PEtabParameter(:ins_app)
Pe_ins_elim =PEtabParameter(:ins_elim)


params =[Pe_a1,Pe_a2,Pe_a3,Pe_a4,Pe_a5]
measurements = DataFrame(obs_id = "obs_xp", time = vec(GlucoDataTime), measurement =vec(val_gluco_dallaMan))
petab_model = PEtabModel(model_mix5comp, obs_xp, measurements, params  ,speciemap  = CI,parametermap = parameter_map)
petab_problem = PEtabODEProblem(petab_model)
p0 = 0.01*ones(length(params))

#res = calibrate(petab_problem,p0,IPNewton())


#print(get_ps(res, petab_problem))

plot(sol_mix.t,glucose_dallaMan(t)/1.88,label="glucose de Dalla-Man",ylims=(0,400),yticks=0:50:400,ylabel="glucose dans le plasma en mg/dL",xlabel = "temps en min ")   #  
plot!(sol_mix.t,sol_mix[:Xp],label="Glucose Bergman + compartiments  ")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   affichage des resultats   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#plot(sol_mix.t,sol_mix[:Xp])
#plot(GlucoDataTime,GlucoDataMean,l=:scatter,label="Gluco data ")
#plot!(sol_mix.t,sol_mix[:Xp],label="Glucose au cours du temps ",ylims=(75,150))
#print(length(params))

#plot(GlucoDataTime,GlucoDataMean,l=:scatter,label="Gluco data ")


#print(get_ps(res, petab_problem))



