#############################################################################
#############################################################################
# ToDo
#
#############################################################################
#############################################################################


import os
import re
import time
import math
import numpy as np
import pandas as pd
import datetime as dtt
import hydroeval as he

def lin_db(x):
    """linear to dB"""
    return 10*np.log10(x)


def db_lin(x):
    """dB to linear"""
    return 10**(x/10)


def timeseries(dates, data):
    """Returns a matrix (list of type(dates,data)) in the format [dates,data]"""
    
    if len(dates)==len(data):
        return [[dates[i],data[i]] for i in range(len(dates))]
    else: raise ValueError(
        f'dates and data must have same first dimension, but have shapes {np.shape(dates)} and {np.shape(data)}')
    

#############################################################################
# Soil Water Balance + Water Cloud model
#############################################################################

def IRR_WCM(PAR, inputs, user_in):
    """Irrigation model and WCM integration.
    
    Based on minimization of KGE between observed and simulated
    $\sigma^0$ values via PSO (pyswarm) optimization.
    The soil water balance model (IRR) produces an estimate of the soil water
    content WW [%] that is used to simulate $\sigma^0$ by a water cloud
    model (WCM).
    
    Inputs
    ----------
    - PAR: initial guess values for parameters to calibrate
        PAR = [A, B, C, D, W_0, W_max, S_fc, S_w, rho_st, Kc]
    - inputs: input quantities for calibration,
        [d, d_sat, P, IRR_obs, EPOT, WWobs, LAI, t_deg, obs]
    - user_in: user-defined options
        irri = user_in: if user_in=True, irrigation is estimated
        and not taken as an input, else the input observed irrigation
        is used in the soil water balance
    
    Return
    -------
    KGE from hydroeval between sigma0 observed and simulated.
    
    """

    # User input
    irri, units = user_in
    
    # Unpack inputs
    A, B, C, D, W_max, WW_fc, WW_w, rho_st, Kc0 = PAR
    t, t_sat, P, IRR_obs, EPOT, WW_obs, WW_sat, veg, angle, sig0_obs = inputs
    
    # Fixed parameters
    # global WW_fc  # = 0.32 # 0.32
    # global WW_w   # = 0.08 # 0.08
    # global rho_st # = 0.4 # /24 # 0.4
    # global Kc0    # = 1 # 0.05 # 1

    W_fc   = WW_fc*W_max # field capacity [mm]
    W_w    = WW_w*W_max # wilting point [mm]
    ET     = np.array([.0]*len(t)) # evapotranspiration
    Ks     = np.array([.0]*len(t)) # water stress coefficient
    rho    = np.array([.0]*len(t)) # depletion fraction
    PS     = np.array([.0]*len(t)) # deep percolation
    W      = np.array([.0]*len(t)) # water content [mm]
    W[0]   = WW_obs[0]*W_max # initial value of sm [mm]    

    if irri==True: IRR = [.0]*len(t) # water content
    else: IRR = IRR_obs
    
    # Build Kc curve
    # Ref. FAO56, Tables 11, 12
    START = t[0].dayofyear; END = t[-1].dayofyear
    START_INI = 134 # first day with vegetation
    L_INI=30; L_IM=40; # length initial and initial to mid stages
    L_MID=45; L_ML=30; # length mid and mid to late stages
    Kc_ini = 0.6*Kc0; Kc_mid = 1.15*Kc0; Kc_late= 0.7*Kc0
    
    # 24 added for compatibility with hourly dataset
    m1 =(Kc_mid-Kc_ini)/L_IM/24 # hourly Kc growth in DEV phase 
    m2 =(Kc_late-Kc_mid)/L_ML/24 # hourly Kc growth in senescence
    
    # Hourly array of Kc values
    Kc_array = np.array([])
    Kc_array = np.append(Kc_array, np.ones((START_INI-START)*24)*Kc_ini)
    Kc_array = np.append(Kc_array, np.ones(L_INI*24)*Kc_ini)
    Kc_array = np.append(Kc_array, [Kc_ini+m1*i for i in range(L_IM*24)])
    Kc_array = np.append(Kc_array, np.ones(L_MID*24)*Kc_mid)
    Kc_array = np.append(Kc_array, [Kc_mid+m2*i for i in range(L_ML*24)])
    Kc_array = np.append(Kc_array, np.ones(len(t)-len(Kc_array))*Kc_late)
    
    for i in [i+1 for i in range(len(t)-1)]:
        
        DOY=t[i].dayofyear
        
        # Build Ks curve
        Kc = Kc_array[i]
        rho[i]=rho_st+0.04*(5-Kc*EPOT[i])
        if W[i-1]>=(1-rho[i])*W_fc:
            Ks[i]=1
        elif (W[i-1]>W_w)and(W[i-1]<(1-rho[i])*W_fc):
            Ks[i]=float(W[i-1]-W_w)/((1-rho[i])*(W_fc-W_w))
        else: Ks[i]=0
        
        # Adjusted evapotranspiration
        ET[i] = EPOT[i]*Kc*Ks[i]
        
        # Irrigation estimate (for summer season only)
        # Irrigation is estimated as the amount of water needed from the day
        # before to take water content up to field capacity
        if irri==True:
            if np.logical_and(DOY>START,DOY<START+100): # summer season
                if W[i-1]<=(1-rho[i])*W_fc: IRR[i]=W_fc-W[i-1]
        
        
        # Water balance [mm]
        W[i]=W[i-1]+P[i]+IRR[i]-ET[i]
        
        # Computation of deep percolation (water above field capacity)
        if W[i]>W_fc:
            PS[i]=W[i]-W_fc
            W[i]=W_fc
            
    WW=np.array(W)/W_max
    WWsat = np.array([ x[1] for x in timeseries(t,WW) if x[0] in t_sat ])
    
    # Water Cloud Model    
    sig0,KGE = WCM([A,B,C,D], [WWsat,veg,angle,sig0_obs], units=units)

    return [WW,IRR,sig0,KGE]


#############################################################################
# Model with all parameters to input by default
#############################################################################

def IRR_WCM_allpar(PAR, inputs, user_in):
    """Irrigation model and WCM integration.
    
    Based on minimization of KGE between observed and simulated
    $\sigma^0$ values via PSO (pyswarm) optimization.
    The soil water balance model (IRR) produces an estimate of the soil water
    content WW [%] that is used to simulate $\sigma^0$ by a water cloud
    model (WCM).
    
    Inputs
    ----------
    - PAR: initial guess values for parameters to calibrate
        PAR = [A, B, C, D, W_0, W_max, S_fc, S_w, rho_st, Kc]
    - inputs: input quantities for calibration,
        [d, d_sat, P, IRR_obs, EPOT, WWobs, LAI, t_deg, obs]
    - user_in: user-defined options
        irri = user_in: if user_in=True, irrigation is estimated
        and not taken as an input, else the input observed irrigation
        is used in the soil water balance
    
    Return
    -------
    KGE from hydroeval between sigma0 observed and simulated.
    
    """

    # User input
    irri, units = user_in
    
    # Unpack inputs
    A, B, C, D, W_max, WW_fc, WW_w, rho_st, Kc0 = PAR
    d, d_sat, P, IRR_obs, EPOT, W_d, W_h, veg, angle, sig0_obs = inputs
    
    # Fixed parameters
    W_0    = W_d[0]
    
    theta  = angle*np.pi/180. # angle of incidence
    W_fc   = WW_fc*W_max # field capacity [mm]
    W_w    = WW_w*W_max # wilting point [mm]
    Ks     = [0]*len(d) # daily, water stress coefficient
    rho    = [0]*len(d) # daily, depletion fraction
    PS     = [0]*len(d) # daily, deep percolation
    W      = [0]*len(d) # daily, water content [mm]
    W[0]   = W_0*W_max # initial value of sm [mm]
    if irri==True: IRR = [0]*len(d) # daily, water content
    else: IRR = IRR_obs
    
    for t in [i+1 for i in range(len(d)-1)]:
        
        DOY=d[t].dayofyear
        
        # Build Kc curve
        # Ref. FAO56, Tables 11, 12
        START  = 134
        Kc_ini = 0.6*Kc0
        Kc_mid = 1.15*Kc0
        Kc_late= 0.7*Kc0

        INI =[x for x in range(START, START+30)] # initial phenological phase
        DEV =[x for x in range(START+30, START+30+40)] # development phase
        m1  =(Kc_mid-Kc_ini)/40 # daily Kc growth in DEV phase 
        MID =[x for x in range(START+30+40, START+30+40+45)] # mid phase
        LATE=[x for x in range(START+30+40+45, START+30+40+45+30)] # late stage
        m2  =(Kc_late-Kc_mid)/30 # daily Kc growth in senescence
        
        if DOY in INI: Kc=Kc_ini
        elif DOY in DEV: Kc=Kc_ini+m1*(DOY-DEV[0])
        elif DOY in MID: Kc=Kc_mid
        elif DOY in LATE: Kc=Kc_mid+m2*(DOY-LATE[0])
        elif DOY>LATE[-1]: Kc=Kc_late
        else: Kc=Kc_ini
        
        # Build Ks curve
        rho[t]=rho_st+0.04*(5-Kc*EPOT[t])
        if W[t-1]>=(1-rho[t])*W_fc:
            Ks[t]=1
        elif (W[t-1]>W_w)and(W[t-1]<(1-rho[t])*W_fc):
            Ks[t]=float(W[t-1]-W_w)/((1-rho[t])*(W_fc-W_w))
        else: Ks[t]=0
        
        # Irrigation estimate (for summer season only)
        # Irrigation is estimated as the amount of water needed from the day
        # before to take water content up to field capacity
        if irri==True:
            if np.logical_and(DOY>START,DOY<START+100): # summer season
                if W[t-1]<=(1-rho[t])*W_fc: IRR[t]=W_fc-W[t-1]
        
        # Water balance [mm]
        W[t]=W[t-1]+P[t]+IRR[t]-EPOT[t]*Kc*Ks[t]
        
        # Computation of deep percolation (water above field capacity)
        if W[t]>W_fc:
            PS[t]=W[t]-W_fc
            W[t]=W_fc
            
    WW=np.array(W)/W_max
    WWsat = np.array([ x[1] for x in timeseries(t,WW) if x[0] in t_sat ])
    
    # Water Cloud Model    
    sig0,KGE = WCM([A,B,C,D], [WWsat,veg,angle,sig0_obs], units=units)

    return [WW,IRR,sig0,KGE]



#############################################################################
# Water Cloud Model #############################################################################

def WCM(PAR, data_in, units='lin'):
    """Water Cloud Model.
    
    This function simulates backscattering with WCM and returns
    the KGE index to perform its minimization for calibration
    of parameters A,B,C,D.
    WCM is parametrized with a single vegetation descriptor (nominated
    LAI, but can be anything).
    Fitting can be performed in linear or dB scale.
    The model is written in one line in both units: this makes the calibration slightly slower, but provides more stability in the parameters' distributions.
    
    
    Inputs
    ------
    - PAR: list
        List of initial guesses for the parameters to calibrate.
    - data_in: list
        List of inputs of observables, that must be in the form:
        [SM,LAI,t_deg,obs], being SM = soil moisture,
        LAI = Leaf Area Index, t_deg = angle of observation,
        obs = observed total sigma0
    - units: str, default 'linear'
        choose to calibrate the model's parameters in 'linear' or 'db' scale
        
    Return
    ------
    KGE between simulated and observed backscattering.
    
    """

    A,B,C,D = PAR # parameters to fit
    WWsat,veg,angle,obs = data_in # input data
    
    theta  = angle*np.pi/180. # angle of incidence
    
    if units=='lin':
        sig0 = lin_db((np.exp((-2*B*veg)/np.cos(theta)))*(db_lin(C+D*WWsat))\
                      +A*veg*np.cos(theta)\
                      *(1-(np.exp((-2*B*veg)/np.cos(theta)))))
    elif units=='db':
        sig0 = (np.exp((-2*B*veg)/np.cos(theta)))*(C+D*WWsat)\
        +(A*veg*np.cos(theta)\
          *(1-(np.exp((-2*B*veg)/np.cos(theta)))))
    else: raise NameError('Please choose one of the options: lin/db')
        
    OUT=he.evaluator(he.kge, sig0, obs) # OUT is kge, r, alpha, beta
    KGE=OUT[0,:][0]

    return [sig0,KGE]

#----------------------------------------------------------------------------

def SM_fromWCM(PAR, data_in, units='lin'):
    """Inverted WCM for SM estimation."""

    A,B,C,D = PAR # parameters, fitted
    WWsat,veg,angle,obs = data_in # input data
    
    theta = angle*np.pi/180. # angle of observation
    T2 = np.exp((-2*B*veg)/np.cos(theta)) # attenuation
    sig0v = A*veg*np.cos(theta)*(1-T2) # sigma0_veg
    
    if units=='lin':
        sig0s_lin = (db_lin(obs)-sig0v)/T2
        SMretr = (lin_db(sig0s_lin)-C)/D
    elif units=='db':
        sig0s = (obs-sig0v)/T2
        SMretr = (sig0s-C)/D 
    
    OUT=he.evaluator(he.kge, SMretr, WWsat) # OUT is kge, r, alpha, beta
    KGE=OUT[0,:];

    return [SMretr,KGE]



#############################################################################
# Soil Water Balance
#############################################################################

def SWB(PAR_SWB, inputs, user_in):
    """Soil water balance + irrigation estimation model .
    
    The soil water balance model (IRR) produces an estimate of the soil water
    content WW [%] and irrigation requirement IRR [mm].
    Ref: Rolle, 2021
    
    Inputs
    ----------
    - PAR: initial guess values for parameters to calibrate
        PAR = [A, B, C, D, W_0, W_max, S_fc, S_w, rho_st, Kc]
    - inputs: input quantities for calibration,
        [d, d_sat, P, IRR_obs, EPOT, WWobs, LAI, t_deg, obs]
    - user_in: user-defined options
        irri = user_in: if user_in=True, irrigation is estimated
        and not taken as an input, else the input observed irrigation
        is used in the soil water balance
    
    Return
    -------
    [WW, IRR]: simulated water content and irrigation requirement 
    
    """

    # User input
    irri, units = user_in
    
    # Unpack inputs
    W_max, WW_fc, WW_w, rho_st, Kc0 = PAR_SWB
    d, P, IRR_obs, EPOT, W_d = inputs
    
    # Fixed parameters
    # global WW_fc  # = 0.32 # 0.32
    # global WW_w   # = 0.08 # 0.08
    # global rho_st # = 0.4 # /24 # 0.4
    # global Kc0    # = 1 # 0.05 # 1
    W_0 = W_d[0]
    
    W_fc   = WW_fc*W_max # field capacity [mm]
    W_w    = WW_w*W_max # wilting point [mm]
    Ks     = [0]*len(d) # daily, water stress coefficient
    rho    = [0]*len(d) # daily, depletion fraction
    PS     = [0]*len(d) # daily, deep percolation
    W      = [0]*len(d) # daily, water content [mm]
    W[0]   = W_0*W_max # initial value of sm [mm]
    if irri==True: IRR = [0]*len(d) # daily, water content
    else: IRR = IRR_obs
    
    for t in [i+1 for i in range(len(d)-1)]:
        
        DOY=d[t].dayofyear
        
        # Build Kc curve
        # Ref. FAO56, Tables 11, 12
        START  = 134 # 14 May
        Kc_ini = 0.6*Kc0
        Kc_mid = 1.15*Kc0
        Kc_late= 0.7*Kc0

        INI =[x for x in range(START, START+30)] # initial phenological phase
        DEV =[x for x in range(START+30, START+30+40)] # development phase
        m1  =(Kc_mid-Kc_ini)/40 # daily Kc growth in DEV phase 
        MID =[x for x in range(START+30+40, START+30+40+45)] # mid phase
        LATE=[x for x in range(START+30+40+45, START+30+40+45+30)] # late stage
        m2  =(Kc_late-Kc_mid)/30 # daily Kc growth in senescence
        
        if DOY in INI: Kc=Kc_ini
        elif DOY in DEV: Kc=Kc_ini+m1*(DOY-DEV[0])
        elif DOY in MID: Kc=Kc_mid
        elif DOY in LATE: Kc=Kc_mid+m2*(DOY-LATE[0])
        elif DOY>LATE[-1]: Kc=Kc_late
        else: Kc=Kc_ini
        
        # Build Ks curve
        rho[t]=rho_st+0.04*(5-Kc*EPOT[t])
        if W[t-1]>=(1-rho[t])*W_fc:
            Ks[t]=1
        elif (W[t-1]>W_w)and(W[t-1]<(1-rho[t])*W_fc):
            Ks[t]=float(W[t-1]-W_w)/((1-rho[t])*(W_fc-W_w))
        else: Ks[t]=0
        
        # Irrigation estimate (for summer season only)
        # Irrigation is estimated as the amount of water needed from the day
        # before to take water content up to field capacity
        if irri==True:
            if np.logical_and(DOY>START,DOY<START+100): # summer season
                if W[t-1]<=(1-rho[t])*W_fc: IRR[t]=W_fc-W[t-1]
        
        # Water balance [mm]
        W[t]=W[t-1]+P[t]+IRR[t]-EPOT[t]*Kc*Ks[t]
        
        # Computation of deep percolation (water above field capacity)
        if W[t]>W_fc:
            PS[t]=W[t]-W_fc
            W[t]=W_fc
            
    WW=np.array(W)/W_max

    return [WW,IRR]



#############################################################################
# Snippets and old versions of the code
#############################################################################

    # Water Cloud Model
    # Multi-line, double units version
    #---------------------------------
    # sig0s_dB = C+D*WWsat # sigma0_soil [dB]
    # T2 = np.exp((-2*B*veg)/np.cos(angle)) # attenuation
    # sig0v = A*veg*np.cos(angle)*(1-T2) # sigma0_veg [units]
    # 
    # if units=='lin':
    #     sig0s = db_lin(sig0s_dB) # sigma0_soil [lin]
    #     sig0_lin = T2*sig0s+sig0v # sigma0_tot [lin]
    #     sig0=lin_db(sig0_lin) # sigma0_tot [dB]
    # elif units=='db':
    #     sig0 = T2*sig0s_dB+sig0v # sigma0_tot [db]
    # else: raise NameError('Please choose one of the options: lin/db')
    
    
    # Old method for WW_sat cut on d_sat dates
    #-----------------------------------------
    # WWsat = pd.DataFrame(timeseries(d,WW))\
    #           .set_index(0).loc[d_sat][1].values