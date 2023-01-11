#!/usr/bin/env python
# coding: utf-8

# # Ausiliary functions (your @staticmethods!)

# Base
import os
import re
import ee
import time
import math
import numpy as np
import pandas as pd
import datetime as dtt

# Analysis
import pyswarms as ps
from scipy import special as sp
from scipy.stats import norm, gamma, f, chi2
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial
from scipy.signal import savgol_filter as sfilter

# Geospatial
import fiona
import xarray as xr
import hydroeval as he
# import geopandas as gpd
# from maps_original import *

# Graphics
import seaborn as sns
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import IPython.display as disp
get_ipython().run_line_magic('matplotlib', 'inline')


# # Statistics and data cleaning / normalization
#-----------------------------------------------------------------------------

def lin_db(x):
    return 10*np.log10(x)

def db_lin(x):
    return 10**(x/10)

#-----------------------------------------------------------------------------

def norm(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

#-----------------------------------------------------------------------------

def bias(obs, sim):
    if len(obs)==len(sim):
        return np.mean(obs-sim)
    else: raise ValueError(
        f'obs and sim must have same first dimension, but have shapes {np.shape(obs)} and {np.shape(sim)}')

#-----------------------------------------------------------------------------

def HIST_norm(ref_mean, ref_std, obs:list):
    """HIST normalization
    Ref. Mladenova, 2013, https://ieeexplore.ieee.org/document/6264094
    
    obs = [value, mean, std]
    """
    value, mean, std = obs
    return ref_mean+ref_std/std*(value-mean)

#-----------------------------------------------------------------------------

def timeseries(dates, data):
    """Returns a matrix (list of type(dates,data)) in the format [dates,data]"""
    
    if len(dates)==len(data):
        return [[dates[i],data[i]] for i in range(len(dates))]
    else: raise ValueError(
        f'dates and data must have same first dimension, but have shapes {np.shape(dates)} and {np.shape(data)}')

#-----------------------------------------------------------------------------


# # Data analysis, fit
# ## Fitting functions

#-----------------------------------------------------------------------------

def linear(x,a,b):
    return a+b*x
    
#-----------------------------------------------------------------------------

def gauss(x, A, mean, dev):
    """Not-normalized, shifted gaussian distribution."""
    
    import math
    
    pdf = (1/(dev*np.sqrt(2*math.pi)))*np.exp(-(x-mean)**2/(2*dev**2))
    return A*pdf

#-----------------------------------------------------------------------------

def skew_gauss(x, A, mean, dev, alpha,):
    """Skew, not-normalized and shifted gaussian distribution.

    References:
    - https://www.wolframalpha.com/input?i=skew+gaussian+distribution
    - https://stackoverflow.com/questions/15400850/scipy-optimize-curve-fit-unable-to-fit-shifted-skewed-gaussian-curve
    - https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewnorm.html

    """
    
    import math
    import scipy.special as sp
    
    pdf = (1/(dev*np.sqrt(2*np.pi)))*np.exp(-pow((x-mean),2)/(2*pow(dev,2)))
    cdf = sp.erfc((-alpha*(x-mean))/(dev*np.sqrt(2)))
    return A*pdf*cdf


#-----------------------------------------------------------------------------
# GEE functions
#-----------------------------------------------------------------------------

def extract_data(image:ee.Image):
    """Ausiliary function to extract data from an Image
    
    This function extracts spatial means and std.dev
    via spatial reducers (reduceRegion).
    Optimal implementation is to map this function
    on a whole ImageCollection via .map() and insert the
    return into a ee.FeatureCollection.
    
    Return
    ------
    ee.Feature
    
    """
    try: # be aware that this try doesn't do anything
        mean = image.reduceRegion(**{ 
            'reducer': ee.Reducer.mean(),
            'geometry': aoi,
        })
        
        dev = image.reduceRegion(**{ 
            'reducer': ee.Reducer.stdDev(),
            'geometry': aoi,
        })
        
        var = image.reduceRegion(**{
            'reducer':ee.Reducer.variance(),
            'geometry': aoi,
        })
    
        properties = {
            'Date': image.get('system:time_start'), # only way to get a timestr is an external operation
            'Geometry': geometry_title,
            'VV[lin]': mean.get('VV'),
            'VH[lin]': mean.get('VH'),
            'Angle[°]': mean.get('angle'),
            'VV_var[lin]': var.get('VV'),
            'VH_var[lin]': var.get('VH'),
            'Orb': image.get('relativeOrbitNumber_start'),
            'Pass': image.get('orbitProperties_pass'),
        }
    except (HttpError, EEException):
        print(f'This image ({image.get("ID")}) had missing data. Skip...\n')
            
            
    return ee.Feature(None, properties)