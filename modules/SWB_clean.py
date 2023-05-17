import sys
sys.path.append('../')

import os
import re
import time
import math
import numpy as np
import pandas as pd
import datetime as dtt
import hydroeval as he
from modules.hallikainen import *


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
# Soil Water Balance
#############################################################################

def SWB(inputs):
    """Soil water balance model
    
    Inputs
    ----------
    - inputs: input quantities for calibration, all as lists or numpy arrays
        WW_fc   = field capacity [m^3/m^3]
        WW_w    = wilting point [m^3/m^3]
        t       = list of dates (datetime object)
        P       = precipitation [mm]
        EPOT    = potential evapotranspiration [mm/day]
        Kc_list = crop coefficient [-], must be computed with function Kc_curve
        rho_st  = standard deplection fraction [-] (set default to 0.5)
        depth   = depth of soil layer under investigation [mm]
    
    Return
    -------
    List of water content (WW), irrigation (IRR), 
    """
    
    # Unpack inputs
    WW_fc, WW_w, t, P, EPOT, Kc_list, rho_st, depth = inputs
    
    Ks      = 0. # water stress coefficient
    rho     = 0. # depletion fraction
    WW      = np.array([.0]*len(t)) # water content [m3/m3]
    WW[0]   = WW_obs[0] # initial value of sm [m3/m3]
    
    for i in [i+1 for i in range(len(t)-1)]:
        
        # Build Ks curve    
        # Compute crop coeff and depletion fraction
        rho=rho_st+0.04*(5-Kc_list[i]*EPOT[i])
        
        if WW[i-1]>=(1-rho)*WW_fc:
            Ks=1
        elif (WW[i-1]>WW_w)and(WW[i-1]<(1-rho)*WW_fc):
            Ks=float(float(WW[i-1]-WW_w)/float((1-rho)*(WW_fc-WW_w)))
        else: Ks=0
        
        if np.logical_and(DOY>121,DOY<244): # summer season between 1 May, 1 Sept
            if W[i-1]<=(1-rho)*W_fc: IRR[i]=W_fc-W[i-1]
        
        # Water balance [mm]
        WW[i]=WW[i-1]+(P[i]+IRR[i]-EPOT[i]*Kc_list[i]*Ks)/(depth)
        
        # Computation of deep percolation (water above field capacity)
        if WW[i]>WW_fc: WW[i]=WW_fc
                
    return [WW,IRR]



#############################################################################
# Kc curves, daily
#############################################################################

# Build Kc curve
# Ref. FAO56, Tables 11, 12

def Kc_curve(t, START_INI, lenghts, coeffs):
    """Ausiliary function to build crop coefficient curves
    
    Follow tables 11, 12 of FAO-56 paper for default values of lenghts of
    phenological stages and Kc values.
    
    Inputs
    ------
    t         = list of dates (datetime object)
    START_INI = first day with vegetation
    lenghts   = list of lenghts of phenological stages [days],
        [initial stage, initial to mid, mid stage, mid to late, late stage]
    coeffs    = list of values of crop coefficient for initial, mid and
        late stage
        
    Returns
    -------
    Crop coefficient curve, as a list of values.
    
    """
    
    L_INI, L_IM, L_MID, L_ML = lenghts
    Kc_ini, Kc_mid, Kc_late = coeffs

    START = t[0].dayofyear; END = t[-1].dayofyear
    m1 =(Kc_mid-Kc_ini)/L_IM # daily Kc growth in DEV phase 
    m2 =(Kc_late-Kc_mid)/L_ML # daily Kc growth in senescence

    # Hourly array of Kc values
    Kc_list = np.array([])
    Kc_list = np.append(Kc_list, np.ones((START_INI-START))*Kc_ini)
    Kc_list = np.append(Kc_list, np.ones(L_INI)*Kc_ini)
    Kc_list = np.append(Kc_list, [Kc_ini+m1*i for i in range(L_IM)])
    Kc_list = np.append(Kc_list, np.ones(L_MID)*Kc_mid)
    Kc_list = np.append(Kc_list, [Kc_mid+m2*i for i in range(L_ML)])
    Kc_list = np.append(Kc_list, np.ones(len(t)-len(Kc_list))*Kc_late)
    
    return Kc_list