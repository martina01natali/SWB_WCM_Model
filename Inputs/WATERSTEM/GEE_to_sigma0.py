#!/usr/bin/env python
# coding: utf-8

# # $\sigma^0$ from GEE + normalization
# 
# This code is a semi-automatic routine for the download of values of $\sigma^0$ from Sentinel-1 via the Google Earth Engine API. The downloaded data represent the mean value of $\sigma^0$ over the AoI, that must be provided manually into the code as a Java polygon or similar object; the variance of the data is provided as well.
# 
# ---
# **Dependencies** 
# 
# This code requires the installation of the Earth Engine API, `ee`. You can find more info on the installation procedure here: [Python installation of GEE](https://developers.google.com/earth-engine/guides/python_install). \
# This code runs on browser-based notebooks only (Google Colaboratory, Jupyter Notebooks, etc...). \
# Be aware that you won't need to install the Google Cloud APK to run the code. 
# 
# ---
# **Normalization**
# 
# After download, $\sigma^0$ values are normalized with respect to orbits. The chosen approach is based on the normalization of the values' distributions of different orbits on a reference distribution (Mladenova, 2013, DOI: 10.1109/TGRS.2012.2205264). The orbit that is chosen as reference by default is the one which average incidence angle is the closest to 40°.
# 
# ---
# **Output**
# 
# The output database has columns:
# - Date: timestamp of date and time of passage rounded at hour
# - Angle[°]: angle of incidence
# - Geometry: name of AoI
# - Orb: orbit relative number
# - Pass: direction of passage, ascending or descending
# - VV[dB], VH[dB]: mean values of $\sigma^0$, not normalized
# - VV_var[dB], VH_var[dB]: variance of values of $\sigma^0$
# - VV_norm[dB], VH_norm[dB]: mean values of $\sigma^0$, normalized

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
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial
from scipy.stats import norm, gamma, f, chi2
from scipy.signal import savgol_filter as sfilter

# Graphics
import seaborn as sns
import matplotlib as mplt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import IPython.display as disp
get_ipython().run_line_magic('matplotlib', 'inline')

# Trigger the authentication flow.
ee.Authenticate()
 
# Initialize the library.
ee.Initialize()


#-----------------------------------------------------------------------------
# Ausiliary functions
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

def gauss(x, A, mean, dev):
    """Not-normalized, shifted gaussian distribution."""
    
    import math
    
    pdf = (1/(dev*np.sqrt(2*math.pi)))*np.exp(-(x-mean)**2/(2*dev**2))
    return A*pdf

#-----------------------------------------------------------------------------

def HIST_norm(ref_mean, ref_std, obs:list):
    """HIST normalization
    Ref. Mladenova, 2013, https://ieeexplore.ieee.org/document/6264094
    
    obs = [value, mean, std]
    """
    value, mean, std = obs
    return ref_mean+ref_std/std*(value-mean)

#-----------------------------------------------------------------------------

# Extract data
print('Mean is computed by spatial average in linear scale.\n'+
      'Std is the square root of variance in linear scale, '+
      'transformed in dB by mantaining constant relative error.')

#-----------------------------------------------------------------------------

def lin_db(x):
    return 10*np.log10(x)

def db_lin(x):
    return 10**(x/10)

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

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Filters definition

# missing data on 2016-10-01 raises not-exceptable exception
# tot = ee.Filter.date('2014-10-03', '2022-12-01')

tot1 = ee.Filter.date('2014-10-03', '2016-09-30')
tot2 = ee.Filter.date('2016-10-02', '2022-12-01')

# Define area of interest
# If you have a GeoJSON file, copy paste.
# If you have a KML, export to GeoJSON (plenty of free tools online)
# or retrieve P

# geoJSON is sanlorenzo2 as an example
geoJSON = {"type":"FeatureCollection", "features": [
{"type":"Feature","geometry":{"type":"Polygon","coordinates":[[[11.347035987013012,43.951592805194856],[11.347393649350675,43.95121129870136],[11.346869077922104,43.95095696103902],[11.346718064935093,43.950774155844215],[11.34676575324678,43.95059929870136],[11.347163155844182,43.95046418181824],[11.347282376623403,43.95048007792214],[11.347139311688338,43.95064698701305],[11.347123415584441,43.95076620779227],[11.347401597402623,43.95084568831175],[11.34753671428574,43.951020545454604],[11.347647987013012,43.95101259740266],[11.347870532467557,43.9508377402598],[11.353148038961045,43.954644857142895],[11.352846012987019,43.95487535064939],[11.352607571428578,43.954692545454584],[11.352416818181826,43.95466075324679],[11.352337337662346,43.954764077922114],[11.352472454545463,43.95495483116887],[11.35244066233767,43.9550263636364],[11.352138636363645,43.95519327272731],[11.351741233766244,43.955368129870166],[11.351502792207803,43.95527275324679],[11.351224610389623,43.955439662337696],[11.35095437662339,43.955415818181855],[11.34934092207794,43.95366724675329],[11.34804538961041,43.952236597402646],[11.347035987013012,43.951592805194856]]]},"properties":{"id":1}}
]}
nfeatures = 1
coords = [geoJSON['features'][i]['geometry']['coordinates'] for i in range(nfeatures)]
aoi = ee.Geometry.MultiPolygon(coords)

geometry_title = input('Please provide a title for AoI geometry. (Default: Budrio_half-right)')
if not geometry_title: geometry_title='Budrio_half-right'


# Get collection of images and filter
img1 = (ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') # linear scale for mean, var computation
        .filterBounds(aoi)
        .filter(tot1)
        .sort('system:time_start'))
img2 = (ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') # linear scale for mean, var computation
        .filterBounds(aoi)
        .filter(tot2)
        .sort('system:time_start'))


def clean_date(date:int):
    return time.strftime('%x %H', time.localtime((date)/1000))

dflist = []
for img in [img1, img2]:
    data = ee.FeatureCollection(img.map(extract_data))
    data_out = data.getInfo()
    data_out_to_df = [e.get('properties') for e in data_out.get('features')]; data_out_to_df[0]
    df = pd.DataFrame.from_dict(data_out_to_df)
    df.Date = df.Date.apply(lambda x : pd.to_datetime(clean_date(x)))
    df['VV[dB]'] = df['VV[lin]'].apply(lambda x : lin_db(x))
    df['VH[dB]'] = df['VH[lin]'].apply(lambda x : lin_db(x))
    df['VV_var[dB]'] = df['VV_var[lin]']/df['VV[lin]']*(10/np.log(10))
    df['VH_var[dB]'] = df['VH_var[lin]']/df['VH[lin]']*(10/np.log(10))
    dflist.append(df)

dftot = pd.concat([dflist[0],dflist[1]]).set_index('Date'); dftot


# # Backscattering normalization
# Ref: DOI: 10.1109/tgrs.2012.2205264, https://ieeexplore.ieee.org/document/6264094 (Mladenova (2013))


import warnings

# need orbits in databaes: exploit sets' uniquity
orb_list = [*set(dftot.Orb.values)]

n=0
leng_list=[]
for orb in orb_list:
    leng = len(dftot[dftot.Orb==orb]) 
    print(f'orbit {orb} has {leng} data')
    if leng==1:
        warnings.warn(f'\nOrbit {orb} has only 1 datum. '+
                      'This orbit will be ignored.')
    leng_list.append(leng)
    n+=leng
print(f'total number of data is {n}')

bad_orb_list=[orb_list[i] for i in range(len(orb_list)) if leng_list[i]==1]
orb_list=[orb_list[i] for i in range(len(orb_list)) if leng_list[i]>1]; orb_list

for bad in bad_orb_list:
    dftot=dftot.drop(dftot[dftot.Orb==bad].index)

normdict = dict()
nnormdict = dict()

for orb in orb_list:
    d = dict()
    for pol in ['VV','VH']:
        d[f'{pol}']      = dftot.loc[dftot.Orb==orb][f'{pol}[dB]'].values
        d[f'{pol}_mean'] = np.mean(dftot.loc[dftot.Orb==orb][f'{pol}[dB]'].values)
        d[f'{pol}_std']  = np.std(dftot.loc[dftot.Orb==orb][f'{pol}[dB]'].values)
        if d[f'{pol}_std']==0: d[f'{pol}_std']=np.sqrt(abs(d[f'{pol}_mean']))
        d[f'angle_mean'] = np.mean(dftot.loc[dftot.Orb==orb][f'Angle[°]'].values)
    d['angle'] = dftot.loc[dftot.Orb==orb][f'Angle[°]'].values
    normdict[orb] = d
    
normdf = pd.DataFrame.from_dict(normdict, orient='index'); normdf
statsnorm = normdf.drop(columns=['VV','VH','angle']); statsnorm

# Find orbit that has mean angle nearest to 40°
# and choose it as reference statistics
orb_ref = normdf.angle_mean.apply(lambda x : abs(x-40)).sort_values().head(1).index.values[0]
if not input(f'The orbit with mean incidence angle nearest to 40° is orbit {orb_ref}. '+
             'This orbit will be taken as reference orbit for normalization. Proceed? [[]/any]')=='':
    raise NameError('Manually set orb_ref variable with chosen reference orbit (type:int).')
else: print(f'\nOrbit {orb_ref} has been set as reference orbit.')


for orb in orb_list:
    d = dict()
    for pol in ['VV', 'VH']:
        d[f'{pol}']      = HIST_norm(
            normdf.loc[orb_ref][f'{pol}_mean'],
            normdf.loc[orb_ref][f'{pol}_std'],
            [
                normdf.loc[orb][pol],
                normdf.loc[orb][f'{pol}_mean'],
                normdf.loc[orb][f'{pol}_std'],
            ]
                 )
        d[f'{pol}_mean'] = np.mean(d[f'{pol}'])
        d[f'{pol}_std']  = np.std(d[f'{pol}'])
    nnormdict[orb] = d
    
nnormdf = pd.DataFrame.from_dict(nnormdict, orient='index'); nnormdf

dftot['VV_norm[dB]']= dftot.apply( 
    lambda x :
    HIST_norm(
        statsnorm['VV_mean'][orb_ref],
        statsnorm['VV_std'][orb_ref],
        [   
            x['VV[dB]'],    
            statsnorm['VV_mean'][x.Orb],    
            statsnorm['VV_std'][x.Orb]
        ],   
    ),
    axis='columns',
)

dftot['VH_norm[dB]']= dftot.apply( 
    lambda x :
    HIST_norm(
        statsnorm['VH_mean'][orb_ref],
        statsnorm['VH_std'][orb_ref],
        [   
            x['VH[dB]'],    
            statsnorm['VH_mean'][x.Orb],    
            statsnorm['VH_std'][x.Orb]
        ],   
    ),
    axis='columns',
)

dftot=dftot.drop(['VH[lin]','VH_var[lin]','VV[lin]','VV_var[lin]'], axis=1)
dftot.drop_duplicates(inplace=True)
dftot.dropna(inplace=True)
dftot

opt_save_df = input('Save df? [y/n]')
if opt_save_df=='y':
    filename = input('Provide filename without extension (def. .csv): [default:geometry name] ')
    if filename=='': filename = geometry_title 
    dftot.to_csv(f'..\Data\{filename}.csv', sep = '\t')

filename = f'..\\Data\\Golden_GEE_2014-22_norm.csv'
if input(f"Wanna save in Data directory? File has name: {filename}. [y/n] ")=='y':
    dftot.to_csv(filename, sep = '\t')

#-----------------------------------------------------------------------------
# ## Plots


def hist_gauss_fit(data, nbins, hist_kwargs, fitline_kwargs,
                   title, density=True, opt_save=False, opt_name='hist_fit',
                  ):
    
    def gauss(x, A, mean, dev):
        """Not-normalized, shifted gaussian distribution."""
        import math
        return A*(1/(dev*np.sqrt(2*math.pi)))*np.exp(-(x-mean)**2/(2*dev**2))

    counts, bins, pads = plt.hist(data, bins=nbins, density=True, **hist_kwargs)
    fit_bounds = [ [0,min(bins),0], [sum(counts)*np.diff(bins)[0],max(bins),abs(max(bins)-min(bins))] ]
    popt, pcov = curve_fit(gauss, bins[:-1], counts, method='trf',bounds=fit_bounds, maxfev=1000)
    A, mean, dev = popt[0], popt[1], popt[2]
    x = np.linspace(min(data), max(data), 50)
    fit = gauss(x, A, mean, dev)
    plt.plot(x, fit, **fitline_kwargs)
    ylabel = 'Density' if density else 'Counts';    plt.ylabel(ylabel)
    plt.legend(loc='best')
    plt.title(title)
    if opt_save: plt.savefig(opt_name+'.png', dpi=300)
    return counts, bins, pads


# Histograms plots

opt_save=True if input('Save plots? [y/n]')=='y' else False

normdf.dropna(inplace=True)
nnormdf.dropna(inplace=True)

# Not normalized
plt.figure()
for orb in orb_list:
    hist_kwargs={'alpha':.5, 'label':f'orb {orb}'}
    fitline_kwargs={'linestyle':'-', 'label':f'fit_{orb}'}
    data = normdf['VV'][orb]
    hist_gauss_fit(data, nbins=10, hist_kwargs=hist_kwargs, fitline_kwargs=fitline_kwargs,
                   title='not normalized', density=True, opt_save=opt_save, opt_name='hist_2014-22_fit_not-norm',)
    plt.xlabel(r'$\sigma^0\ [VV]$ [dB]')
t = plt.text(np.min(bins)+1, 0.9*np.max(counts),
             f'# entries: {len(dftot.index)}',
             ha="center", va="center", size=15,
             bbox=dict(boxstyle="round,pad=0.3", fc="tab:orange", ec="k", lw=2, alpha=.5))

# Normalized
plt.figure()
for orb in orb_list:
    hist_kwargs={'alpha':.5, 'label':f'orb {orb}'}
    fitline_kwargs={'linestyle':'-', 'label':f'fit_{orb}'}
    data = nnormdf['VV'][orb]
    counts, bins, pads = hist_gauss_fit(data, nbins=10, hist_kwargs=hist_kwargs, fitline_kwargs=fitline_kwargs,
                   title='normalized', density=True, opt_save=opt_save, opt_name='hist_2014-22_fit_norm',)
    plt.xlabel(r'$\sigma^0\ [VV]$ [dB]')
t = plt.text(np.min(bins)+1, 0.9*np.max(counts),
             f'# entries: {len(dftot.index)}',
             ha="center", va="center", size=15,
             bbox=dict(boxstyle="round,pad=0.3", fc="tab:orange", ec="k", lw=2, alpha=.5))


# Timeseries, before normalization

plot_name = 'sigma_norm_2014-22_'
opt_save = True if input('Save plot? [y/n]')=='y' else False
plt.figure(figsize=(20,5))

for orb in orb_list:
    plt.plot(dftot[dftot.Orb==orb]['VV[dB]'], alpha=.5, marker='o', linestyle='-', label=f'orb {orb}')
plt.ylabel(r'$\sigma^0$ [dB]')
plt.xlabel('index')
plt.legend(loc='best')
plt.title('not normalized')
if opt_save: plt.savefig(plot_name+'plot_not-norm.png', dpi=300)
plt.show()


# Timeseries, after normalization

plot_name = 'sigma_norm_2014-22_'
opt_save = True if input('Save plot? [y/n]')=='y' else False

plt.figure(figsize=(20,5))
for orb in orb_list:
    plt.plot(dftot[dftot.Orb==orb].index, dftot[dftot.Orb==orb]['VV[dB]'], alpha=.5, marker='o', linestyle='-', label=f'orb {orb}')
plt.ylabel(r'$\sigma^0$ [dB]')
plt.xlabel('index')
plt.legend(loc='best')
plt.title('normalized')
if opt_save: plt.savefig(plot_name+'plot_norm.png', dpi=300)
plt.show()