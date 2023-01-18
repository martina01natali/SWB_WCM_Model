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
    

#----------------------------------------------------------------------------

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
    A, B, C, D, W_max, Kc = PAR
    d, d_sat, P, IRR_obs, EPOT, W_d, W_h, veg, theta, VV = inputs
    
    # Fixed parameters
    W_0 = W_d[0]
    WW_fc = 0.32
    WW_w = 0.08
    rho_st = 0.4    
    
    angle = theta*np.pi/180. # angle of incidence
    W_fc = WW_fc*W_max # field capacity [mm]
    W_w = WW_w*W_max # wilting point [mm]
    Ks = [0]*len(d) # daily, water stress coefficient
    rho = [0]*len(d) # daily, depletion fraction
    PS = [0]*len(d) # daily, deep percolation
    W = [0]*len(d) # daily, water content [mm]
    W[0] = W_0*W_max # initial value of sm [mm]
    if irri==True: IRR = [0]*len(d) # daily, water content
    else: IRR = IRR_obs
    
    for t in [i+1 for i in range(len(d)-1)]:
        rho[t]=rho_st+0.04*(5-Kc*EPOT[t])
        if W[t-1]>=(1-rho[t])*W_fc:
            Ks[t]=1
        elif (W[t-1]>W_w)and(W[t-1]<(1-rho[t])*W_fc):
            Ks[t]=float(W[t-1]-W_w)/((1-rho[t])*(W_fc-W_w))
        else: Ks[t]=0
        
        DOY=d[t].dayofyear
        
        # Irrigation estimate (for summer season only)
        # Irrigation is estimated as the amount of water needed from the day
        # before to take water content up to field capacity
        if irri==True:
            if np.logical_and(DOY>134,DOY<235): # summer season
                if W[t-1]<=(1-rho[t])*W_fc: IRR[t]=W_fc-W[t-1]
        
        # Water balance [mm]
        W[t]=W[t-1]+P[t]+IRR[t]-EPOT[t]*Kc*Ks[t]
        
        # Computation of deep percolation (water above field capacity)
        if W[t]>W_fc:
            PS[t]=W[t]-W_fc
            W[t]=W_fc
            
    WW=np.array(W)/W_max
    WWsat = pd.DataFrame(timeseries(d,WW)).set_index(0).loc[d_sat][1].values
    
    
    sig0s_dB = C+D*WWsat # sigma0_soil [dB]
    T2 = np.exp((-2*B*veg)/np.cos(angle)) # attenuation
    
    if units=='lin':
        sig0s = db_lin(sig0s_dB) # sigma0_soil [lin]
        sig0v = A*veg*np.cos(angle)*(1-T2) # sigma0_veg [lin]
        sig0_lin = T2*sig0s+sig0v # sigma0_tot [lin]
        sig0=lin_db(sig0_lin) # sigma0_tot [dB]
    elif units=='db':
        sig0v = A*veg*np.cos(angle)*(1-T2) # sigma0_veg [db]
        sig0 = T2*sig0s+sig0v # sigma0_tot [db]
    else: raise NameError('Please choose one of the options: lin/db')

    
    # T2 = np.exp((-2*B*veg)/np.cos(angle)) # two-way attenuation from the vegetation layer
    # sig0s = db_lin(C+D*WWsat) # define bare soil backscatter [fit in dB, then in lin]
    # sig0v = A*veg*np.cos(angle)*(1-T2) # define backscatter from the vegetation [fit in lin]
    # sig0_lin = T2*sig0s+sig0v
    # sig0=lin_db(sig0_lin) # from linear scale to dB
        
    OUT=he.evaluator(he.kge, sig0, VV) # OUT is kge, r, alpha, beta
    KGE=OUT[0,:];

    return [WW,IRR,sig0,KGE]