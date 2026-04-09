function Sb(pm_5,pm_6,pHE_b) # from eq (6)
    (pm_6-pHE_b)/pm_5
end

function m_4(pm_5,pm_6,pIb,pV_I,pHE_b) # from eq (9)
    (2/5)*(Sb(pm_5,pm_6,pHE_b)/(pIb*pV_I))*(1-pHE_b)
end

function m_2(pm_5,pm_6,pIb,pV_I,pHE_b) # from eq (9)
    
    (Sb(pm_5,pm_6,pHE_b)/(pIb*pV_I)-(m_4(pm_5,pm_6,pIb,pV_I,pHE_b)/(1-pHE_b)))*(1-pHE_b)/pHE_b
end

function Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2) # from eq (20)
    (pF_cns-pEGPb+pk_1*pGb*pV_G)/pk_2
end
    
function V_m0(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2,pxK_m0) # from eq (22)
     (pEGPb-pF_cns)*(pxK_m0+Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2))/Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2)
 end

function k_p1(pk_p2,pk_p3,pk_p4,pEGPb,pGb,pV_G,pγ,pSb,pIb) # from eq (12)
    pEGPb + pk_p2*pGb*pV_G+pk_p3*pIb+pk_p4*pSb/pγ
end

I(I_p)  = I_p/V_I                         # plasma insulin concentration
HE(S)  = -m_5 * S + m_6                   # hepatic extraction
m_3(S)  = (HE(S) * m_1)/(1.0-HE(S))       # hepatic degradation

m_3b(HE_b,m_1)    = (HE_b*m_1)/(1.0-HE_b) # basal hepatic degradation

I_bf() = I_b

# algebraic functions
Ra_meal(Q_gut) = (f*k_abs*Q_gut)/BW  # exogenous glucose rate of appearance, mg/kg/min
Ex_meal(Q_gut) = k_abs*Q_gut*(1.0-(f/BW))
k_empt(Q_sto,Dose) = k_min + ((k_max-k_min)/2.0) * (tanh(α(Dose)*(Q_sto - b*Dose)) - tanh(β(Dose)*(Q_sto - d * Dose)) + 2.0) # gastric emptying, /min
α(Dose) = 5.0/(2.0*Dose*(1.0-b)) # rate of decrease to k_min, /mg
β(Dose) = 5.0/(2.0*Dose*d)       # rate of increase to k_max, /mg
@syms Q_gut(t)
Ra_meal_time(t) = Ra_meal(Q_gut(t))


G_emptEla(t,Dose) = (Dose*β_ela*(k^β_ela)*(t^(β_ela-1.0))*exp(-(k*t)^β_ela))
Ra_ela(Q_gut) = (f*k_abs*Q_gut)/BW
Ex_ela(Q_gut) = k_abs*Q_gut *(1.0-(f/BW))

egp(G_p,id,I_po) = k_p1(k_p2,k_p3,k_p4,EGP_b,G_b,V_G,γ,Sb(m_5,m_6,HE_b),I_b) - k_p2*G_p - k_p3*id - k_p4* I_po # endogenous glucose production, mg/kg/k_min


U_ii()      = F_cns  # insulin independent utilization                                                        
U_id(X,G_t) = V_m(X) / (K_m(X) + G_t) # insulin-dependent glucose utilization, mg/kg/min

V_m(X)      = V_m0(F_cns,EGP_b,k_1,G_b,V_G,k_2,K_m0) + V_mx*X
K_m(X)      = K_mx*X + K_m0



E(G_p) = ifelse(G_p>k_e2,k_e1*(G_p - k_e2),0.0)

heaviside(x,func) = func(0.5*(sign(x)+1.0)) # return 0 or 1 , x: value , func : floor or ceil function
Ysec(Y,G_p) = IfElse.ifelse(β_Y*((G_p/V_G)-h)-(-Sb(m_5,m_6,HE_b))>=0,(-α_Y*Y)+(((α_Y*β_Y)*((G_p/V_G)-h))), (-α_Y*Y)-(α_Y*Sb(m_5,m_6,HE_b)))
z(Y,G_p,G_t,id,I_po,Q_gut) = IfElse.ifelse(dG(G_p,G_t,id,I_po,Q_gut) > 0.0, Y+K*dG(G_p,G_t,id,I_po,Q_gut)+Sb(m_5,m_6,HE_b), Y+Sb(m_5,m_6,HE_b))
dG(G_p,G_t,id,I_po,Q_gut) = (egp(G_p,id,I_po) + Ra_meal(Q_gut) - U_ii() - E(G_p) - (k_1 * G_p) + (k_2 * G_t))/V_G 
