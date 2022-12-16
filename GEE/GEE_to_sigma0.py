#!/usr/bin/env python
# coding: utf-8

import ee
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm, gamma, f, chi2
import IPython.display as disp
get_ipython().run_line_magic('matplotlib', 'inline')

# Trigger the authentication flow.
ee.Authenticate()
 
# Initialize the library.
ee.Initialize()

#-----------------------------------------------------------------------------
# Define area of interest
# If you have a GeoJSON file, copy paste.
# If you have a KML, export to GeoJSON (plenty of free tools online)
# or retrieve P

geoJSON = {"type":"FeatureCollection", "features": [
{"type":"Feature","geometry":{"type":"Polygon","coordinates":[[[11.347035987013012,43.951592805194856],[11.347393649350675,43.95121129870136],[11.346869077922104,43.95095696103902],[11.346718064935093,43.950774155844215],[11.34676575324678,43.95059929870136],[11.347163155844182,43.95046418181824],[11.347282376623403,43.95048007792214],[11.347139311688338,43.95064698701305],[11.347123415584441,43.95076620779227],[11.347401597402623,43.95084568831175],[11.34753671428574,43.951020545454604],[11.347647987013012,43.95101259740266],[11.347870532467557,43.9508377402598],[11.353148038961045,43.954644857142895],[11.352846012987019,43.95487535064939],[11.352607571428578,43.954692545454584],[11.352416818181826,43.95466075324679],[11.352337337662346,43.954764077922114],[11.352472454545463,43.95495483116887],[11.35244066233767,43.9550263636364],[11.352138636363645,43.95519327272731],[11.351741233766244,43.955368129870166],[11.351502792207803,43.95527275324679],[11.351224610389623,43.955439662337696],[11.35095437662339,43.955415818181855],[11.34934092207794,43.95366724675329],[11.34804538961041,43.952236597402646],[11.347035987013012,43.951592805194856]]]},"properties":{"id":1}}
]}

nfeatures = 1 # number of polygons to include
coords = [geoJSON['features'][i]['geometry']['coordinates'] for i in range(nfeatures)]
aoi = ee.Geometry.MultiPolygon(coords)

geometry_title = input('Please provide a title for AoI geometry. (Default: Budrio_half-right)')
if not geometry_title: geometry_title='Budrio_half-right'

#-----------------------------------------------------------------------------
# Filters definition

sp17 = ee.Filter.date('2017-04-04', '2017-05-22') # Budrio-related filters
su17 = ee.Filter.date('2017-05-22', '2017-09-15') # Budrio-related filters
au17 = ee.Filter.date('2017-09-15', '2017-11-02') # Budrio-related filters
tot = ee.Filter.date('2014-10-03', '2022-12-01')

# Missing data on 2016-10-01 raises not-exceptable exception
# Solution: two date filters, concatenate databases
tot1 = ee.Filter.date('2014-10-03', '2016-09-30')
tot2 = ee.Filter.date('2016-10-02', '2022-12-01')

#-----------------------------------------------------------------------------
# Get collection of images and filter

img1 = (ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') # linear scale for mean, var computation
        .filterBounds(aoi)
        .filter(tot1)
        .sort('system:time_start'))
img2 = (ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') # linear scale for mean, var computation
        .filterBounds(aoi)
        .filter(tot2)
        .sort('system:time_start'))

#-----------------------------------------------------------------------------
# Extract data

print('Mean is computed by spatial average in linear scale.\n'+
      'Std is the square root of variance in linear scale, '+
      'transformed in dB by mantaining constant relative error.')


def lin_db(x):
    return 10*np.log10(x)

def db_lin(x):
    return 10**(x/10)

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

dftot = pd.concat([dflist[0],dflist[1]]).set_index('Date')

opt_save_df = input('Save df? [y/n]')
if opt_save_df=='y':
    filename = input('Provide filename without extension (def. .csv): ')
    dftot.to_csv(f'..\Data\{filename}.csv', sep = '\t')