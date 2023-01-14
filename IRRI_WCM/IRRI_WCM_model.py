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
    backscattering via PSO (pyswarm) optimization.
    Version with 8 parameters to calibrate.
    
    Inputs
    ----------
    - PAR: initial guess values for parameters to calibrate
        PAR = [A, B, C, D, W_0, W_max, S_fc, S_w, rho_st, Kc]
    - inputs: input quantities for calibration,
        [d, d_sat, P, IRRobs, EPOT, WWobs, LAI, t_deg, obs]
    
    Return
    -------
    KGE from hydroeval between sigma0 observed and simulated.
    
    """

    # User input
    irri = user_in
    
    # Unpack inputs
    A, B, C, D, W_0, W_max, S_fc, S_w, rho_st, Kc = PAR
    d, d_sat, P, IRRobs, EPOT, WWobs, LAI, t_deg, obs = inputs
    
    #S_fc = .46 # hardcoded
    #S_w = .08 # hardcoded
    W_fc = S_fc*W_max # water content at field capacity
    W_w  = S_w*W_max # water content at wilting point
    theta = t_deg*np.pi/180. # angle of observation
    
    
    if irri==True: IRR = [0]*len(d) # daily, water content
    else: IRR = IRRobs
    
    Ks = [0]*len(d) # daily, water stress coefficient
    rho = [0]*len(d) # daily, depletion fraction
    PS = [0]*len(d) # daily, deep percolation
    W = [0]*len(d) # daily, water content
    
    W[0] = W_0*W_max
    
    for t in [i+1 for i in range(len(d)-1)]:
        rho[t]=rho_st+0.04*(5-Kc*EPOT[t])
        if W[t-1]>=(1-rho[t])*W_fc:
            Ks[t]=1
        elif (W[t-1]>W_w)and(W[t-1]<(1-rho[t])*W_fc):
            Ks[t]=float(W[t-1]-W_w)/((1-rho[t])*(W_fc-W_w))
        else: Ks[t]=0
        
        DOY=d[t].dayofyear
        
        # Irrigation estimate (for summer season only)
        # Irrigation is estimated as the amount of water needed at the day before
        # to take water content up to field capacity
        if np.logical_and(DOY>134,DOY<235): # summer season
            if irri==True:
                if W[t-1]<=(1-rho[t])*W_fc: IRR[t]=W_fc-W[t-1]
        
        # Water balance
        W[t]=W[t-1]+P[t]+IRR[t]-EPOT[t]*Kc*Ks[t]
        
        # Computation of deep percolation (water above field capacity)
        if W[t]>W_fc:
            PS[t]=W[t]-W_fc
            W[t]=W_fc
            
    WW=np.array(W)/W_max   
    WWsat = pd.DataFrame(timeseries(d,WW)).set_index(0).loc[d_sat][1].values
    
    T2 = np.exp((-2*B*LAI)/np.cos(theta)) # two-way attenuation from the vegetation layer
    sig0s = db_lin(C+D*WWsat) # define bare soil backscatter [fit in dB, then in lin]
    sig0v = A*LAI*np.cos(theta)*(1-T2) # define backscatter from the vegetation [fit in lin]
    sig0_lin = T2*sig0s+sig0v
    sig0=lin_db(sig0_lin) # from linear scale to dB
        
    OUT=he.evaluator(he.kge, sig0, obs) # OUT is kge, r, alpha, beta
    KGE=OUT[0,:];

    return [WW,IRR,sig0,KGE]