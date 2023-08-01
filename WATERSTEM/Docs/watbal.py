def singleLayerSWB(w1, int_pio, et0):
    """Simple water balance

    f_min: minimum water content
    coeff_p: related to capillarity
    k_sat: saturated hydraulic conductivity
    coeff_c: related to soil compaction

    sat_eff: saturation efficiency
    int_pio: potential infiltration
    ric: hydraulic conductivity
    AET: actual evapotranspiration
    spe_str: specific storage, or soil moisture storage capacity
    pas_tem_sec: time step
    """

    global w_sat, w_fc, w_w, f_min, coeff_p, k_sat, coeff_c, spe_str, pas_tem_sec

    sat_eff = w1/w_sat # saturation efficiency
    cap_inf = f_min+coeff_p*(1-sat_eff) # capillary infiltration
    inf = int_pio if int_pio <= cap_inf else cap_inf # net infiltration
    ric = k_sat*sat_eff**coeff_c
    stress = w_fc # stress threshold is set to field capacity
    
    if w1>=w_fc:
        AET=et0
    elif (w1<w_fc)&(w1>w_w):
        AET=et0*(w1-w_w)/(stress-w_w)
    else:
        AET=0

    w2=w1+( (inf - ric - AET)/spe_str*pas_tem_sec )
    w2=w2 if w2<=w_sat else w_sat # cut at w_sat for runoff 
    
    return w2


def 2layerSWB(PIO, EPOT, SWE, W, W2, Ks, m2, Ks2, m22, W_max, W_max2, alpha):
    """soil water balance over 2 layers

    """

    IE = (PIO + SWE) * ((W / W_max) ** alpha)

    E1 = EPOT * W / W_max
    E2 = 0
    E = E1 + E2

    PERC = Ks * (W / W_max) ** m2 * (W2 < W_max2)
    PERC_amm = W + PIO + SWE - IE - E1
    PERC_amm = PERC_amm * (PERC_amm > 0)
    PERC = PERC * (PERC <= PERC_amm) + PERC_amm * (PERC > PERC_amm)
    W = max(0, W + PIO + SWE - IE - PERC - E1)

    PERC2 = Ks2 * (W2 / W_max2) ** m22
    PERC2_amm = W2 + PERC - E2
    PERC2_amm = PERC2_amm * (PERC2_amm > 0)
    PERC2 = PERC2 * (PERC2 <= PERC2_amm) + PERC2_amm * (PERC2 > PERC2_amm)
    W2 = max(0, W2 + PERC - PERC2 - E2)

    SE = (W - W_max) * (W > W_max)
    SE2 = (W2 - W_max2) * (W2 > W_max2)
    W[W > W_max] = W_max
    W2[W2 > W_max2] = W_max2
    W[W < 0] = 0
    W2[W2 < 0] = 0

    BF = PERC2
    QS = IE + SE + SE2

    return BF, QS, W, W2, E



