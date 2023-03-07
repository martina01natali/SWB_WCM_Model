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
from scipy import stats
from scipy.stats import norm, gamma, f, chi2
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial
from scipy.signal import savgol_filter as sfilter

# Geospatial
import hydroeval as he

# Graphics
import seaborn as sns
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates



#############################################################################
# Fitting functions
#############################################################################

def gauss(x, sigmag, mu, A):
    """Not-normalized, shifted gaussian distribution."""
        
    pdf = (1/(sigmag*np.sqrt(2*math.pi)))*np.exp(-(x-mu)**2/(2*sigmag**2))
    return A*pdf


def skew_gauss(x, sigmag, mu, alpha, a, c,):
    
    normpdf = (1/(sigmag*np.sqrt(2*math.pi)))*np.exp(-(np.power((x-mu),2)/(2*np.power(sigmag,2))))
    normcdf = (0.5*(1+sp.erf((alpha*((x-mu)/sigmag))/(np.sqrt(2)))))
    
    return 2*a*normpdf*normcdf + c



#############################################################################
# Histogram with fit
#############################################################################

def hist_gauss_fit(data, nbins, hist_kwargs, fitline_kwargs,
                   title, density=False,
                   opt_save=False, dir_name='', opt_name='hist_fit',
                   func='optim', thr_asymm=2, fit_method='trf',
                  ):
    """Histogram with automatic gaussian fit.
    
    Arguments
    ---------
    - func: object, default gauss
        WARNING: skew_gauss not supported yet
        
    """

    mean = np.mean(data); rang=np.ptp(data)
    x = np.linspace(mean-rang, mean+rang, 1000)
    counts, bins, pads = plt.hist(data, bins=nbins, density=density, **hist_kwargs)
    
    thr_asymm = thr_asymm
    if func=='optim':
        if abs(stats.skew(data))<thr_asymm: func_fit=gauss
        else: func_fit=skew_gauss
    else: func_fit=func
        
    integ = sum(np.diff(bins)*counts)/sum(counts) if density==False else 1
    
    if func_fit==gauss:
        fit_label='gauss'
        kwargs = {
        'bounds': [ 
            [0,min(data),0],
            [np.std(data)*2,max(data),integ*2] ]
        }
    elif func_fit==skew_gauss:
        fit_label='skew_norm'
        kwargs = {
            'p0' : (np.std(data), np.mean(data),
                    stats.skew(data), integ, min(data)),
            'bounds' : [
                [0, min(0,np.mean(data)*2),  min(0,stats.skew(data)*2), min(0,integ*2), min(0,min(data))],
                [np.std(data)*2, max(0,np.mean(data)*2), max(0,stats.skew(data)*2), max(0,integ*2), max(0,min(data))]
            ]
        }
    else: raise ValueError(f'Func {func} is not a valid option.')
    
    
    popt, pcov = curve_fit(func_fit, bins[:-1], counts, method=fit_method, **kwargs, maxfev=100000)
    fit = func_fit(x, *popt)
    plt.plot(x, fit, label=fit_label+' fit', **fitline_kwargs)
    ylabel = 'Density' if density else 'Counts'; plt.ylabel(ylabel)
    plt.xlabel(title);
    plt.legend(loc='best');
    
    if opt_save: plt.savefig(dir_name+opt_name+'.png', dpi=300)
    
    return [counts, bins, pads, popt, pcov]


#############################################################################
# Analysis on parameters out of multiple runs
#############################################################################

# def parameters_analysis(params, PAR, PARn_str, bounds, timestr, automate, nbins=10):
#     
#     matrix = np.array(
#         [
#             np.array(
#                 [ params[i][j] for i in range(len(params)) ])
#             for j in range(len(PAR))
#         ]
#     )
#     
#     PAR_dict = {
#         'A':      ['[-]',        round(A, 2) ],
#         'B':      ['[-]',        round(B, 2) ],
#         'C':      ['[dB]',       round(C, 2) ],
#         'D':      ['[dB m3/m3]', round(D, 2) ],
#         'W_max':  ['[mm]',       round(W_max, 2) ],
#         'WW_fc':  ['[m3/m3]',    round(WW_fc, 2) ],
#         'WW_w':   ['[m3/m3]',    round(WW_w, 2) ],
#         'rho_st': ['[mm/h]',     round(rho_st, 2) ],
#         'Kc0' :   ['[-]',        round(Kc0, 2)],
#     }
# 
#     PARn = []
#     PARn_dev = []
#     rowLabels = []
#     
#     if automate: opt_save=True
#     else: opt_save = True if input('Save histograms of params? [y/n]')=='y' else False
#     
#     for label in PAR_dict:
#         if label in PARn_str:
#             i = PARn_str.index(label)
#             data = matrix[i]
#             
#             if len(data)==1:
#                 PARn=[x[0] for x in matrix]
#                 PARn_dev=[0 for x in matrix]
#             else:
#                 hist_kwargs={'alpha':.5, }
#                 fitline_kwargs={'linestyle':'-',}
#                 counts, bins, pads, popt, pcov = hist_gauss_fit(
#                     data, nbins=nbins, hist_kwargs=hist_kwargs,
#                     fitline_kwargs=fitline_kwargs,
#                     title=f'{label} {PAR_dict[label][0]}', density=True,
#                     opt_save=opt_save, dir_name='Plot\\',
#                     opt_name=f'{timestr}_hist_{label}',
#                     func='optim', thr_asymm=2, fit_method='dogbox')
#                 mean=np.mean(data); rang=np.ptp(data)
#                 plt.xlim(mean-rang, mean+rang);
#                 
#                 PARn.append(popt[1])
#                 PARn_dev.append(popt[0])
#                 plt.show()
#             
#             PAR_dict[label].pop(1)
#             PAR_dict[label].append('cal');
#             PAR_dict[label].append([round(bounds[0][i],2), round(bounds[1][i],2)]);
#             PAR_dict[label].append(np.round(PARn[i], 3));
#             PAR_dict[label].append(np.round(PARn_dev[i], 3));
#             PAR_dict[label].append(matrix[i]);            
#             
#             if not 'Cal/fix' in rowLabels: rowLabels.append('Cal/fix')
#             if not 'Bounds' in rowLabels: rowLabels.append('Bounds')
#             if not 'Mean' in rowLabels: rowLabels.append('Mean')
#             if not 'Dev' in rowLabels: rowLabels.append('Dev')
#             
#         else:
#             PAR_dict[label].append('fix')
#             PAR_dict[label].append(['/', '/'])
#             PAR_dict[label].append(PAR_dict[label][1])
#             PAR_dict[label].append('/')
#             PAR_dict[label].pop(1)
#     
#     print(PAR_dict, '\n', rowLabels)
#     return PAR_dict, rowLabels