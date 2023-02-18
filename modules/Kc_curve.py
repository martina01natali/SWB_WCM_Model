# Kc curves

import os
import re
import time
import math
import numpy as np
import pandas as pd
import datetime as dtt


# Build Kc curve
# Ref. FAO56, Tables 11, 12

def Kc_curve(year, t):

    if year==2017:
        START = t[0].dayofyear; END = t[-1].dayofyear
        START_INI = pd.Timestamp(2017, 5, 23).dayofyear-30 # first day with vegetation
        L_INI=30; L_IM=40; # length initial and initial to mid stages
        L_MID=45; L_ML=30; # length mid and mid to late stages
        Kc_ini = 0.6; Kc_mid = 1.15; Kc_late= 0.7
    
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
        
    elif year==2020:
        START = t[0].dayofyear; END = t[-1].dayofyear
        START_INI = pd.Timestamp(2020, 3, 25).dayofyear # first day with vegetation
        ADD=30 # days to add to take into account initial growth until anthesis
        L_INI=ADD+20; L_IM=25; # length initial and initial to mid stages
        L_MID=25; L_ML=10; # length mid and mid to late stages
        Kc_ini = 0.3; Kc_mid = 1.15; Kc_late= 1.05
    
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
    
    return Kc_array