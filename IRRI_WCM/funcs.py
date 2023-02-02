#!/usr/bin/env python
# coding: utf-8

# Ausiliary functions (your @staticmethods!)
# To import all the functions in this module, run:
#
# from dir.file import *
#
# being dir the directory of file.


# Base
import os
import re
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

#----------------------------------------------------------------------------# Statistics and data cleaning / normalization

def lin_db(x):
    """linear to dB"""
    return 10*np.log10(x)


def db_lin(x):
    """dB to linear"""
    return 10**(x/10)


def norm(x):
    """min-max normalization of x data"""
    return (x-np.min(x))/(np.max(x)-np.min(x))


def norm_fit(x, a, b):
    """norm in range a,b
    
    x = data, a<b
    """
    return a+(x-np.min(x))*(b-a)/(np.max(x)-np.min(x))


def bias(obs, sim):
    """distance between obs' and sim's mean values"""
    if len(obs)==len(sim):
        return np.mean(obs-sim)
    else: raise ValueError(
        f'obs and sim must have same first dimension, but have shapes {np.shape(obs)} and {np.shape(sim)}')

    
def Rvalue(x:list,y:list)->float:
    """compute Pearson's R between x,y data"""
    if len(x)!=len(y):
        raise ValueError(
            f'x and y must have same first dimension,'+
            'but have shapes{np.shape(x)} and {np.shape(y)}')
        
    matrix = np.array( 
        [ [x[i], y[i]] for i in range(len(x))
         if not (np.isnan(x[i])and(np.isnan(y[i]))) ] )
    return np.corrcoef(matrix,rowvar=False)[0][1]

    
def HIST_norm(ref_mean, ref_std, obs:list):
    """HIST normalization
    Ref. Mladenova, 2013, https://ieeexplore.ieee.org/document/6264094
    
    obs = [value, mean, std]
    """
    value, mean, std = obs
    return ref_mean+ref_std/std*(value-mean)

#----------------------------------------------------------------------------
def timeseries(dates, data):
    """Returns a matrix (list of type(dates,data)) in the format [dates,data]"""
    
    if len(dates)==len(data):
        return [[dates[i],data[i]] for i in range(len(dates))]
    else: raise ValueError(
        f'dates and data must have same first dimension, but have shapes {np.shape(dates)} and {np.shape(data)}')
    

#----------------------------------------------------------------------------
# Data analysis, fit
# Fitting functions

def gauss(x, A, mean, dev):
    """Not-normalized, shifted gaussian distribution."""
    import math
    return A*(1/(dev*np.sqrt(2*math.pi)))*np.exp(-(x-mean)**2/(2*dev**2))

def skew_gauss(x, A, mean, dev, alpha):
    """Skew, not-normalized and shifted gaussian distribution."""

    import math
    import scipy.special as sp
    
    pdf = (1/(dev*np.sqrt(2*np.pi)))*np.exp(-pow((x-mean),2)/(2*pow(dev,2)))
    cdf = sp.erfc((-alpha*(x-mean))/(dev*np.sqrt(2)))
    return A*pdf*cdf

def hist_gauss_fit(data, nbins, hist_kwargs, fitline_kwargs,
                   title, density=False,
                   opt_save=False, dir_name='', opt_name='hist_fit',
                   func=gauss,
                  ):
    """Histogram with automatic gaussian fit.
    
    Arguments
    ---------
    - func: object, default gauss
        WARNING: skew_gauss not supported yet    
        
    """

    x = np.linspace(min(data), max(data), 200)
    counts, bins, pads = plt.hist(data, bins=nbins, density=density, **hist_kwargs)
    if func==gauss:
        fit_bounds = [ [0,min(bins),0],
                      [sum(counts)*np.diff(bins)[0],max(bins),abs(max(bins)-min(bins))] ]
    elif func==skew_gauss:
        fit_bounds = [ [0,min(bins),0],
                      [sum(counts)*np.diff(bins)[0],max(bins),abs(max(bins)-min(bins))],
                      # [-10, 10]
                     ]
    else: raise ValueError(f'Func {func} is not a valid option.')
    popt, pcov = curve_fit(func, bins[:-1], counts, method='trf',bounds=fit_bounds, maxfev=1000)
    fit = func(x, *popt)
    plt.plot(x, fit, **fitline_kwargs)
    ylabel = 'Density' if density else 'Counts'; plt.ylabel(ylabel)
    plt.legend(loc='best'); plt.title(title)
    xtext= 0.5*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] # 0.8*(max(x)-min(x))+min(x);
    ytext=0.5*max(counts) #0.5*(max(counts)-min(counts))+min(counts)
    t = plt.text(xtext, ytext,
                 f'tot counts={len(data)}\nmean={popt[1]:.2f}\ndev={popt[2]:.2f} ({popt[2]/abs(popt[1])*100:.1f}%)',
                 ha="center", va="center", size=15,
                 bbox=dict(boxstyle="round,pad=0.3", fc="tab:orange", ec="k", lw=2, alpha=.5))
    
    if opt_save: plt.savefig(dir_name+opt_name+'.png', dpi=300)
    
    return [counts, bins, pads, popt, pcov]


#----------------------------------------------------------------------------
# Ausiliary, file/dir management

