# packages
using Catalyst
using IfElse
using Latexify
using OrdinaryDiffEq
using Plots


@variables t 
@species  G_p(t) G_t(t) I_l(t) Q_sto1(t) Q_sto2(t)  I_1(t)  X(t) Y(t) I_d(t) I_po(t) I_p(t) Q_gut(t)


###################################
# computation of some parameters  #
###################################

function Sb(pm_5,pm_6,pHE_b) # from eq (6)
    (pm_6-pHE_b)/pm_5
end

function m_4(pm_5,pm_6,pIb,pV_I,pHE_b) # from eq (9)
    (2/5)*(Sb(pm_5,pm_6,pHE_b)/(pIb*pV_I))*(1-pHE_b)
end

function m_2(pm_5,pm_6,pIb,pV_I,pHE_b) # from eq (9)
    
    (Sb(pm_5,pm_6,pHE_b)/(pIb*pV_I)-(m_4(pm_5,pm_6,pIb,pV_I,pHE_b)/(1-pHE_b)))*(1-pHE_b)/pHE_b
end

# to compute the basal glucose Gb, one should solve the following equations:
# using Reduce
# eq1 = :(pEGPb = Uidb + pF_cns) # from eq (2) and assuming Eb = 0 (ie. fasting glucose is always lower than excression threshold)
# eq2 = :(Gtb * pV_m0 = (pEGPb-pF_cns)*(pxK_m0 + Gtb)) # from eq (22)
# eq3 = :(pk_1 * Gpb = pEGPb - pF_cns + pk_2*Gtb) # from eq (20)
 #sol = Algebra.solve((eq1,eq2,eq3),(:G_pb,:Uidb,:G_tb))

# function Gb(pk_1, pk_2, pV_m0, pxK_m0, pF_cns, pEGPb, pV_G)
#     ((pk_2 * pxK_m0 + pV_m0 + 2pF_cns) * pEGPb - ((pk_2 * pxK_m0 + pV_m0 + pF_cns) * pF_cns + pEGPb ^ 2)) / (((pF_cns + pV_m0) - pEGPb) * pk_1)/pV_G
# end

# function Gtb(pF_cns,pEGPb,pk_1,pV_G,pk_2,pV_m0,pxK_m0) # from eq (20)
#     pGb = Gb(pk_1, pk_2, pV_m0, pxK_m0, pF_cns, pEGPb)
#     (pF_cns-pEGPb+pk_1*pGb*pV_G)/pk_2
# end

function Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2) # from eq (20)
    (pF_cns-pEGPb+pk_1*pGb*pV_G)/pk_2
end
    
function V_m0(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2,pxK_m0) # from eq (22)
     (pEGPb-pF_cns)*(pxK_m0+Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2))/Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2)
 end

function k_p1(pk_p2,pk_p3,pk_p4,pEGPb,pGb,pV_G,pγ,pSb,pIb) # from eq (12)
    pEGPb + pk_p2*pGb*pV_G+pk_p3*pIb+pk_p4*pSb/pγ
end

##############################################
# Glucose Kinatics: fluxes plasma <-> tissues
##############################################

# parameters 
@parameters V_G k_1 k_2
# algebraic function
G(G_p) = G_p/V_G # plasma glucose concentration, mg/dL

# chemical reactions network
GK = @network_component GK begin
    k_1, G_p --> G_t # glucose: plasma --> tissue, mg/kg --> mg/kg 
    k_2, G_t --> G_p # glucose: tissue --> plasma, mg/kg --> mg/kg
end # k_1 k_2

# GKeq = convert(ODESystem,GK)
# latexify(GKeq) |> render

##############################################
# Insulin Kinatics: fluxes plasma <-> liver
##############################################
@parameters m_1 m_5 m_6 V_I γ I_b HE_b

# algebraic functions
I(I_p)  = I_p/V_I                         # plasma insulin concentration
HE(S)  = -m_5 * S + m_6                   # hepatic extraction
m_3(S)  = (HE(S) * m_1)/(1.0-HE(S))       # hepatic degradation

m_3b(HE_b,m_1)    = (HE_b*m_1)/(1.0-HE_b) # basal hepatic degradation

I_bf() = I_b

# chemical reactions network
IK = @network_component IK begin
    @species I_po(t) 
    m_1,                           I_l  --> I_p  # insulin going from liver to plasma, pmol/kg --> pmol/kg
    m_2(m_5,m_6,I_b,V_I,HE_b),     I_p  --> I_l  # insulin going from plasma to liver, pmol/kg --> pmol/kg
    m_3(γ*I_po),                   I_l  --> ∅    # hepatic degradation of insulin 
    m_4(m_5,m_6,I_b,V_I,HE_b), I_p  --> ∅    # spontaneous degradation of insulin
end # m_1 m_5 m_6 V_I HE_b γ I_b

#IKeq = convert(ODESystem,IK)
#latexify(IKeq) |> render

##### REMARQUE 
##### paramètre V_I inutile ???

#########################################
# Gastro-intestinal tract (Cobelli 2006)
#########################################

@parameters k_gri k_min k_max k_abs f b d BW Dose

# algebraic functions
Ra_meal(Q_gut) = (f*k_abs*Q_gut)/BW  # exogenous glucose rate of appearance, mg/kg/min
Ex_meal(Q_gut) = k_abs*Q_gut*(1.0-(f/BW))
k_empt(Q_sto,Dose) = k_min + ((k_max-k_min)/2.0) * (tanh(α(Dose)*(Q_sto - b*Dose)) - tanh(β(Dose)*(Q_sto - d * Dose)) + 2.0) # gastric emptying, /min
α(Dose) = 5.0/(2.0*Dose*(1.0-b)) # rate of decrease to k_min, /mg
β(Dose) = 5.0/(2.0*Dose*d)       # rate of increase to k_max, /mg
Ra_meal_time(t) = Ra_meal(Q_gut(t))


# chemical reactions network
RA = @network_component RA begin
    @species Q_gut(t)
    k_gri,                 Q_sto1 --> Q_sto2 # glucose from solid to liquid in the stomach, mg --> mg
    k_empt(Q_sto1+Q_sto2,Dose), Q_sto2 --> Q_gut  # glucose going from the stomach to the gut, mg --> mg 
    Ra_meal(Q_gut),        Q_gut ⇒ G_p        # glucose intestinal absorption, mg --> mg/kg
    Ex_meal(Q_gut),        Q_gut ⇒ ∅
end # k_gri k_min k_max k_abs f b d BW Dose

#RAeq = convert(ODESystem,RA)
#latexify(RAeq) |> render

#########################################
# Gastro-intestinal tract (Elashoff 1982)
#########################################

@parameters Dose β_ela k f k_abs BW

G_emptEla(t,Dose) = (Dose*β_ela*(k^β_ela)*(t^(β_ela-1.0))*exp(-(k*t)^β_ela))
Ra_ela(Q_gut) = (f*k_abs*Q_gut)/BW
Ex_ela(Q_gut) = k_abs*Q_gut *(1.0-(f/BW))

ELA = @network_component ELA begin
    G_emptEla(t,Dose), Q_sto1 ⇒ Q_gut
    Ra_ela(Q_gut), Q_gut ⇒ G_p
    Ex_ela(Q_gut), Q_gut ⇒ ∅
end #Dose β_ela k f k_abs BW

#ELAeq = convert(ODESystem,ELA)
#latexify(ELAeq) |> render

#########################################
# Endogenous Glucose Production (EGP)
#########################################

# parameters
@parameters k_p2 k_p3 k_p4 k_i EGP_b G_b I_b

# algebraic functions
egp(G_p,id,I_po) = k_p1(k_p2,k_p3,k_p4,EGP_b,G_b,V_G,γ,Sb(m_5,m_6,HE_b),I_b) - k_p2*G_p - k_p3*id - k_p4* I_po # endogenous glucose production, mg/kg/k_min

# chemical reactions network
EGP = @network_component EGP begin
    @species I_p(t)
    egp(G_p,I_d,I_p), ∅   --> G_p # endogenous glucose production
    k_i,               I_1 --> I_d    # delayed insulin action, pmol/L --> pmol/L
    k_i,               I_d --> ∅      # delayed insulin action, pmol/L --> ∅
    k_i*I_p/V_I,       ∅   --> I_1    # delayed insulin action, pmol/L --> pmol/L
end #k_p2 k_p3 k_p4 k_i V_I m_5 m_6 HE_b

#EGPeq = convert(ODESystem,EGP)
#latexify(EGPeq) |> render

#####################################################
# Glucose Uptake (insulino dependent and independent)
#####################################################

# parameters
@parameters V_mx K_mx K_m0 F_cns p2u I_b

# algebraic functions
U_ii()      = F_cns  # insulin independent utilization                                                        
U_id(X,G_t) = V_m(X) / (K_m(X) + G_t) # insulin-dependent glucose utilization, mg/kg/min

V_m(X)      = V_m0(F_cns,EGP_b,k_1,G_b,V_G,k_2,K_m0) + V_mx*X
K_m(X)      = K_mx*X + K_m0

# chemical reactions network
GUptake = @network_component GUptake begin
    @species I_p(t)
    U_ii(), G_p ⇒ ∅ # insulino-independent glucose utilization
    U_id(X,G_t), G_t --> ∅
    p2u, X --> ∅ 
    p2u*(I_p/V_I), ∅ --> X
    I_b*p2u, X ⇒ ∅
end # p2u V_mx K_m0 K_mx V_I F_cns I_b

#GUptakeEq = convert(ODESystem,GUptake)
#latexify(GUptakeEq) |> render
##################
# Renal excretion
##################

# parameters
@parameters k_e1 k_e2

# algebraic fu"nctions
#(G_p) = ifelse(G_p>k_e2,k_e1*(G_p - k_e2),0.0)
E(G_p) = ifelse(G_p>k_e2,k_e1*(G_p - k_e2),0.0)

# chemical reactions network
RE = @network_component RE begin
   E(G_p), G_p --> ∅ # renal glucose clearance 
end # k_e1 k_e2

#REeq = convert(ODESystem,RE)
#latexify(REeq) |> render

####################
# Insulin Secretion
####################
heaviside(x,func) = func(0.5*(sign(x)+1.0)) # return 0 or 1 , x: value , func : floor or ceil function

# parameters
@parameters γ K α_Y β_Y k_1 k_2

# algebraic functions
####### REMARQUE
h = G_b
## Attention: il manquait un β_Y dans le heavyside:
#Ysec(Y,G_p) = ( heaviside(β_Y*((G_p/V_G)-h)-(-Sb(m_5,m_6,HE_b)),ceil) * ((-α_Y*Y)+(((α_Y*β_Y)*((G_p/V_G)-h)))) ) + ( heaviside((-Sb(m_5,m_6,HE_b))-β_Y*((G_p/V_G)-h),floor) * ((-α_Y*Y)-(α_Y*Sb(m_5,m_6,HE_b))) )
## Let use ifelse instead of heavyside:
Ysec(Y,G_p) = IfElse.ifelse(β_Y*((G_p/V_G)-h)-(-Sb(m_5,m_6,HE_b))>=0,(-α_Y*Y)+(((α_Y*β_Y)*((G_p/V_G)-h))), (-α_Y*Y)-(α_Y*Sb(m_5,m_6,HE_b)))
S_po(Y,G_p,G_t,id,I_po,Q_gut) = IfElse.ifelse(dG(G_p,G_t,id,I_po,Q_gut) > 0.0, Y+K*dG(G_p,G_t,id,I_po,Q_gut)+Sb(m_5,m_6,HE_b), Y+Sb(m_5,m_6,HE_b))
dG(G_p,G_t,id,I_po,Q_gut) = (egp(G_p,id,I_po) + Ra_meal(Q_gut) - U_ii() - E(G_p) - (k_1 * G_p) + (k_2 * G_t))/V_G 
# chemical reactions network
IS = @network_component IS begin
    @species I_d(t) G_p(t) Q_gut(t) G_t(t)
    Ysec(Y,G_p), ∅ --> Y      # insulin secretion in pancreas
    S_po(Y,G_p,G_t,I_d,I_po,Q_gut), ∅ --> I_po   # insulin from pancreas to portal vein  
    γ ,                             I_po --> I_l
 
end # γ K α_Y β_Y k_1 k_2 m_5 m_6 HE_b G_b

#ISeq = convert(ODESystem,IS)
#latexify(ISeq) |> render

#################################
# Insulin Secretion (simplified)
#################################
Ysec2(Y,G_p) = (-α_Y*(Y-β_Y*((G_p/V_G)-h)))

# chemical reactions network
IS2 = @network_component TS2 begin
    @species I_d(t) G_t(t) Q_gut(t) G_p(t)
    Ysec2(Y,G_p),                    ∅ --> Y      # insulin secretion in pancreas
    S_po(Y,G_p,G_t,I_d,I_po,Q_gut), ∅ --> I_po   # insulin from pancreas to portal vein  
    γ ,                             I_po --> I_l
end #γ K α_Y β_Y k_1 k_2 m_5 m_6 HE_b G_b

##########
# DXylose
##########

# parameters
@parameters k_min k_max k_abs f b d BW DoseDX k_1 k_2 DXk_e1 DXk_e2

# algebraic functions
EDX(DX_p) = IfElse.ifelse(DX_p>DXk_e2,(DXk_e1*(DX_p - DXk_e2)^2.0),0.0)

# chemical reaction networks
DX = @network_component DX begin
    k_max, DX_sto1 --> DX_sto2 # gastric grinding
    k_empt(DX_sto1 + DX_sto2,DoseDX), DX_sto2 --> DX_gut # gastric emptying
    k_abs, DX_gut  --> ∅ # intestinal absorption
    Ra_meal(DX_gut), ∅ --> DX_p # intestinal absorption
    k_1, DX_p --> DX_t # blood --> tissue
    k_2, DX_t --> DX_p # tissue --> blood
    EDX(DX_p),DX_p ⇒ ∅ # renal clearance
end # k_min k_max k_abs f b d BW DoseDX k_1 k_2 DXk_e1 DXk_e2 

#DXeq = convert(ODESystem,DX)
#latexify(DXeq) |> render

##############
# DXylose ELA
##############

# parameters
@parameters β_ela k f k_abs BW k_1 k_2 DXk_e1 DXk_e2 DoseDX

DXELA = @network_component DXELA begin
    G_emptEla(t,DoseDX), DX_sto1 ⇒ DX_gut
    Ra_ela(DX_gut), DX_gut ⇒ DX_p
    Ex_ela(DX_gut), DX_gut ⇒ ∅
    
    k_1, DX_p --> DX_t # blood --> tissue
    k_2, DX_t --> DX_p # tissue --> blood
    EDX(DX_p),DX_p ⇒ ∅ # renal clearance    
end # β_ela k f k_abs BW k_1 k_2 DXk_e1 DXk_e2 DoseDX

#DXELAeq =  convert(ODESystem,DXELA)
#latexify(DXELAeq) |> render

#############
# C-peptides
#############

# parameters
@parameters γ k_01 k_21 k_12 V_C BW 

# chemical reaction networks
CP = @network_component CP begin
    ((γ*I_po)/V_C)*BW, ∅ --> CP_1
    k_01, CP_1 --> ∅
    k_21, CP_1 --> CP_2
    k_12, CP_2 --> CP_1
end # γ k_01 k_21 k_12 V_C BW

#CPeq = convert(ODESystem,CP)
#latexify(CPeq) |> render 

include("params.jl")

@named s_m0 = extend(IK , GK )
@named s_m1 = extend(s_m0,RA)
@named s_m2 = extend(s_m1, EGP)
@named s_m3 = extend(s_m2,GUptake)
@named s_m4 = extend(s_m3,RE)
@named full_model = extend(s_m4,IS2)
full_model_complete = complete(full_model)
#full_ode = convert(ODESystem, full_model_complete)

println(full_model)
p = all_p[:diab]
typeof(p)
CI=gen_init_states(p)
CI2 = gen_init_states(all_p[:norm])

print(CI)
print(CI2)
#println(p)

tspan = (0.0, 420.0)
print(CI)
prob_dallaman = ODEProblem(complete(full_model), CI, tspan, p)


sol_dallaman = solve(prob_dallaman)
plot(sol_dallaman)
#plot(sol_dallaman.t,sol_dallaman[:Q_sto2])
#plot(sol_dallaman.t,Ra_meal_time)

plot(sol_dallaman.t,sol_dallaman[:I_p]/0.04,ylims=(0,400),yticks=0:50:400)   
plot!(sol_dallaman.t,sol_dallaman[:G_p]/1.49,ylims=(0,400),yticks=0:50:400)  

