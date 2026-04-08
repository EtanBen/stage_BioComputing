# ODE model definition

function dalla_net!(du, u, p, t)
    # parameters
    pV_G, pk_1, pk_2, pV_I, pm_1, pm_5, pm_6, pHE_b, pk_max, pk_min, pk_abs, pk_gri, pf, pb, pc, pBW, pD, pk_p2, pk_p3, pk_p4, pk_i, pF_cns, pGb, pV_mx, pxK_m0, pxK_mx, pp_2U, pK, pxα_Y, pxβ_Y, pγ, pk_e1, pk_e2, pIb, pEGPb = p

    pSb = Sb(pm_5,pm_6,pHE_b)
    # pGtb = Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2)
    pV_m0 = V_m0(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2,pxK_m0)
    #pGb = Gb(pk_1, pk_2, pV_m0, pxK_m0, pF_cns, pEGPb, pV_G)
    pk_p1 = k_p1(pk_p2,pk_p3,pk_p4,pEGPb,pGb,pV_G,pγ,pSb,pIb)
    pm_4 = m_4(pm_5,pm_6,pIb,pV_I,pHE_b)
    pm_2 = m_2(pm_5,pm_6,pIb,pV_I,pHE_b)
    pα = 5 / (2*pD*(1-pb))
    pβ = 5 / (2*pc*pD)
    ph = pGb

    # algebraic equations
    G     = t -> u[1] / pV_G # G

    I     = t -> u[4] / pV_I                  # I
    m_3   = t -> (HE(t) * pm_1) / (1 - HE(t)) # m_3
    HE    = t -> - pm_5 * S(t) + pm_6          # HE

    Ra     = t     -> (pf * pk_abs * u[7]) / pBW # Ra
    q_sto  = t     -> u[5] + u[6]             # q_sto
    k_empt = q_sto -> pk_min + (( pk_max -  pk_min) / 2.) * (tanh( pα * (q_sto -  pb *  pD)) - tanh( pβ * (q_sto - pc *  pD))+ 2.) # k_empt 

    EGP   = t -> max(pk_p1 - pk_p2 * u[1] - pk_p3 * u[9] - pk_p4 * u[12],0) # EGP

    U_ii   = t   -> pF_cns                                     # U_ii
    U_id   = t   -> (V_m(u[10]) * u[2]) / (K_m(u[10]) + u[2]) # U_id
    V_m    = X_t -> pV_m0 + pV_mx * X_t                         # V_m
    K_m    = X_t -> pxK_m0 + pxK_mx * X_t                         # K_m

    function piecewiseDY(t)          # d.Y
        if pxβ_Y * (G(t) - ph) >= - pSb
            return - pxα_Y * (u[11] - pxβ_Y * (G(t) - ph))
        elseif pxβ_Y * (G(t) - ph) < - pSb
            return - pxα_Y * u[11] - pxα_Y * pSb
        end
    end
    S      = t -> pγ * u[12]          # S
    S_po   = t ->  piecewiseS_po(t)  # S_po
    function piecewiseS_po(t)        # S_po   
        if (du[1] / pV_G) > 0.
            return u[11] + pK * (du[1] / pV_G) + pSb
        elseif (du[1]/ pV_G) <= 0.
            return u[11] + pSb
        end
    end
    # ---------------------
    # Renal Excretion
    E = t -> piecewiseE(t) # E
    function piecewiseE(t) # E
        if u[1] > pk_e2
            return pk_e1 * (u[1] - pk_e2)
        elseif u[1] <= pk_e2
            return 0.
        end
    end


    # ------------------------------------------------------------------------
    # differential equations

    # Glucose kinetics
    du[1] = EGP(t) + Ra(t) - U_ii(t) - E(t) - pk_1 * u[1] + pk_2 * u[2] # d.G_p
    du[2] = - U_id(t) + pk_1 * u[1] - pk_2 * u[2]                       # d.G_t


    # --------------------
    # Insulin kinetics
    du[3] = - (pm_1 + m_3(t)) * u[3] + pm_2 * u[4] + S(t) # d.I_l
    du[4] = - (pm_2 + pm_4) * u[4] + pm_1 * u[3]           # d.I_p


    # -------------------- 
    # Glucose rate of appearance
    #du[5]  = - pk_gri * u[5] + pD * DiracDelta(0.001)(t) # d.q_sto1
    du[5]  = - pk_gri * u[5] # d.q_sto1
    du[6]  = - k_empt(q_sto(t)) * u[6] + pk_gri * u[5]  # d.q_sto2
    du[7]  = - pk_abs * u[7] + k_empt(q_sto(t)) * u[6]  # d.q_gut


    # --------------------
    # Endogenous production
    du[8] = - pk_i * (u[8] - I(t)) # d.I_1
    du[9] = - pk_i * (u[9] - u[8]) # d.I_d


    # --------------------
    # Utilization
    du[10] = - pp_2U * u[10] + pp_2U * (I(t)-pIb) # d.X


    # --------------------
    # Secretion
    du[11] = piecewiseDY(t)          # d.Y
    du[12] = - pγ * u[12] + S_po(t)   # d.I_po    
end # dalla_net!
