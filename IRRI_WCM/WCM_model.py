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

#----------------------------------------------------------------------------

def WCM(PAR, data_in, units='lin'):
    """Water Cloud Model.
    
    This function simulates backscattering with WCM and returns
    the KGE index to perform its minimization for calibration
    of parameters A,B,C,D.
    WCM is parametrized with a single vegetation descriptor (nominated
    LAI, but can be anything).
    Fitting can be performed in linear or dB scale.
    
    
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
    SM,LAI,t_deg,obs = data_in # input data
    
    theta = t_deg*np.pi/180. # angle of observation
    sig0s_dB = C+D*SM # sigma0_soil [dB]
    T2 = np.exp((-2*B*LAI)/np.cos(theta)) # attenuation
    
    if units=='lin':
        sig0s = db_lin(sig0s_dB) # sigma0_soil [lin]
        sig0v = A*LAI*np.cos(theta)*(1-T2) # sigma0_veg [lin]
        sig0_lin = T2*sig0s+sig0v # sigma0_tot [lin]
        sig0=lin_db(sig0_lin) # sigma0_tot [dB]
    elif units=='db':
        sig0v = A*LAI*np.cos(theta)*(1-T2) # sigma0_veg [db]
        sig0 = T2*sig0s+sig0v # sigma0_tot [db]
    else: raise NameError('Please choose one of the options: linear/db')
        
    OUT=he.evaluator(he.kge, sig0, obs) # OUT is kge, r, alpha, beta
    KGE=OUT[0,:];

    return [sig0,KGE]

#----------------------------------------------------------------------------

def SM_fromWCM(PAR, data_in, units='lin'):
    """Inverted WCM for SM estimation."""

    A,B,C,D = PAR # parameters, fitted
    SM,LAI,t_deg,obs = data_in # input data
    
    theta = t_deg*np.pi/180. # angle of observation
    T2 = np.exp((-2*B*LAI)/np.cos(theta)) # attenuation
    sig0v = A*LAI*np.cos(theta)*(1-T2) # sigma0_veg
    
    if units=='lin':
        sig0s_lin = (db_lin(obs)-sig0v)/T2
        SMretr = (lin_db(sig0s_lin)-C)/D
    elif units=='db':
        sig0s = (obs-sig0v)/T2
        SMretr = (sig0s-C)/D 
    
    OUT=he.evaluator(he.kge, SMretr, SM) # OUT is kge, r, alpha, beta
    KGE=OUT[0,:];

    return [SMretr,KGE]