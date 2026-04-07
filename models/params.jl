# parameter values (both normal and diabetic)
DM_params = Dict(
              V_G   => [1.88,1.49], 
              k_1   => [0.065,0.042],  
              k_2   => [0.079,0.071], 
              V_I   => [0.05,0.04], 
              m_1   => [0.190,0.379], 
              # :m_2   => [0.484,0.637], # can be computed from equation (9)
              # :m_4   => [0.194,0.269], # can be computed from equation (9)
              m_5   => [0.0304,0.0526], 
              m_6   => [0.6471,0.8118],
              HE_b  => [0.6,0.6],
              k_max => [0.0558,0.0465], 
              k_min => [0.0080,0.0076], 
              k_abs => [0.057,0.023], 
              k_gri => [0.0558,0.0465],  
              f     => [0.9,0.9], 
              b     => [0.82,0.68], 
              d     => [0.010,0.09], 
              # α     => [0.00013,0.00006], # can be computed from eq (10) in [DallaMan & al. 2006]
              # β     => [0.00236,0.00023], # can be computed from eq (11) in [DallaMan & al. 2006]
              BW    => [78.0,91.0],
              Dose  => [90000.0, 90000.0], 
              # k_p1  => [2.70,3.09], # can be computed from equation (12)
              k_p2  => [0.0021,0.0007], 
              k_p3  => [0.009,0.005], 
              k_p4  => [0.0618,0.0786], 
              k_i   => [0.0079,0.0066], 
              F_cns => [1.0,1.0],
              # V_m0  => [2.50,4.65], # We keep this parameter even if it can be computed as in eq 22 of Dalla Man's paper
              V_mx  => [0.047,0.034], 
              K_m0 => [225.59,466.21], 
              K_mx => [0.0,0.0],
              p2u  => [0.0331,0.0840], 
              K     => [2.30,0.99], 
              α_Y  => [0.050,0.013], 
              β_Y  => [0.11,0.05], 
              γ     => [0.5,0.5], 
              # h     => [89.39980053,141.32659636], # h is supposed = to Gb
              k_e1  => [0.0005,0.0007], 
              k_e2  => [339.0,269.0], 
              I_b    => [25.49, 54.81], 
              # Sb    => [1.54,3.57], # I remove this parameter because it is computed from HEb pSb = (m6-HEb)/m5
              # we add 2 parameters to allow for estimation of basal Glucose and EGP
              G_b    => [91.76, 164.18],
             
              
              EGP_b  => [1.92, 2.01])

all_pn = [V_G,k_1, k_2, V_I, m_1, m_5, m_6, HE_b, k_max, k_min, k_abs, k_gri, f, b, d, BW, Dose, k_p2, k_p3, k_p4, k_i, F_cns, G_b, V_mx, K_m0, K_mx, p2u, K, α_Y, β_Y, γ, k_e1, k_e2, I_b, EGP_b]

# parameters are more convenient to handle with the following dictionary
all_p = Dict(st => Dict(pn => DM_params[pn][i] for pn in all_pn) 
             for (st,i) in zip([:norm,:diab],[1,2]))

# all variable names
all_vn = [G_p,G_t,I_l,I_p,Q_sto1,Q_sto2,Q_gut,I_1,I_d,X,Y,I_po]

"""
    traces_of_obs_vars(sol,                         # solution of an ODE
                       rn::ReactionSystem,          # reaction network
                       all_p::Dict{Num,Float64}     # all parameter values
                      )::Dict{Symbol,Vector{Float64}}  #

return the tuple of the traces of G, I, EGP, Ra, U and S.
"""
function traces_of_obs_vars(sol,rn::ReactionSystem,all_p::Dict{Num,Float64})::Dict{Symbol,Vector{Float64}}
    # pGtb = Gtb(all_p[:F_cns],all_p[:EGP_b],all_p[:k_1],all_p[:G_b],all_p[:V_G],all_p[:k_2])
    pSb = Sb(all_p[m_5],all_p[m_6],all_p[HE_b])
    #pGb = Gb(all_p[k_1],all_p[k_2], all_p[V_m0], all_p[K_m0], all_p[F_cns], all_p[EGP_b], all_p[V_G])
    pGb = all_p[G_b]
    pk_p1 = k_p1(all_p[k_p2],all_p[k_p3],all_p[k_p4],all_p[EGP_b],pGb,all_p[V_G],all_p[γ],pSb,all_p[I_b])
    pV_m0 = V_m0(all_p[F_cns],all_p[EGP_b],all_p[k_1],pGb,all_p[V_G],all_p[k_2],all_p[K_m0])

    Vm = X -> pV_m0 + all_p[V_mx] * X
    Km = X -> all_p[K_m0] + all_p[K_mx] * X

    idxGp = indexin(G_p,species(rn))[1]
    idxIp = indexin(I_p,species(rn))[1]
    idxId = indexin(I_d,species(rn))[1]
    idxIpo = indexin(I_po,species(rn))[1]
    idxQgut = indexin(Q_gut,species(rn))[1]
    idxX = indexin(X,species(rn))[1]
    idxGt = indexin(G_t,species(rn))[1]

    traces = Dict()

    traces[:G] = sol[idxGp,:] / all_p[V_G] 
    traces[:I] = sol[idxIp,:] / all_p[V_I] 
    traces[:EGP] = -all_p[k_p2]*sol[idxGp,:] - all_p[k_p3]*sol[idxId,:] - all_p[k_p4]*sol[idxIpo,:] .+ pk_p1
    # EGP = max.(-all_p[:k_p2]*sol[1,:] - all_p[:k_p3]*sol[9,:] - all_p[:k_p4]*sol[12,:] .+ pk_p1,0.0)
    traces[:Ra] = all_p[f] * all_p[k_abs] * sol[idxQgut,:] / all_p[BW]
    traces[:U] = (Vm.(sol[idxX,:]) .* sol[idxGt,:])./(Km.(sol[idxX,:]) .+ sol[idxGt,:]) .+ all_p[F_cns]
    traces[:S] = all_p[γ] * sol[idxIpo,:]

    traces
end # traces_of_obs_vars

"""
    gen_init_states(params::Dict{Num,Float64})::Dict{Num,Float64}

returns the init states of all variables as a dictionary
"""
function gen_init_states(params::Dict{Num,Float64})::Dict{Num,Float64}
    init_states::Dict{Num,Float64} = Dict()
    m30::Float64 = params[HE_b]*params[m_1]/(1-params[HE_b])
    pSb = Sb(params[m_5],params[m_6],params[HE_b])
    pm_4 = m_4(params[m_5],params[m_6],params[I_b],params[V_I],params[HE_b])
    pm_2 = m_2(params[m_5],params[m_6],params[I_b],params[V_I],params[HE_b])
    # pGb = Gb(params[k_1],params[k_2], params[V_m0], params[K_m0], params[F_cns], params[EGP_b], params[V_G])
    pGb = params[G_b]
    #pV_m0 = pV_m0(all_p[F_cns],all_p[EGP_b],all_p[k_1],pGb,all_p[V_G],all_p[k_2],all_p[K_m0])

    # pSb = (params[m_6]-params[HE_b])/params[m_5]
    # pm_4 = (2/5)*(pSb/(params[I_b]*params[V_I]))*(1-params[HE_b]) # from eq (9)
    # pm_2 = (pSb/(params[I_b]*params[V_I])-(pm_4/(1-params[HE_b])))*(1-params[HE_b])/params[HE_b] # from eq (9)

    init_states[G_p] = pGb*params[V_G]
    init_states[G_t] = (params[F_cns]-params[EGP_b]+params[k_1]*init_states[G_p])/params[k_2]
    init_states[I_p] = params[I_b]*params[V_I]
    init_states[I_l] = init_states[I_p]*(pm_4+pm_2)/params[m_1]
    init_states[Q_sto1] = params[Dose]
    init_states[Q_sto2] = 0.0
    init_states[Q_gut] = 0.0
    init_states[I_1] = params[I_b]
    init_states[I_d] = params[I_b]
    init_states[X] = 0.0
    init_states[Y] = 0.0
    #init_states["Ipo"] = (init_states["Ip"]*pm_4+init_states["Il"]*m30)/params[γ]
    init_states[I_po] = pSb/params[γ] # see eq (23) in Dalla Man's paper

    init_states
end # gen_init_state