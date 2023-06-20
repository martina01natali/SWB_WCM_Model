# Ref. https://pyeto.readthedocs.io/en/latest/hargreaves.html

# ET0 calculation

import datetime
import sys

sys.path.append('../')
from modules.pyeto.pyeto import *

def hargre(lat_deg, dates, temp_min, temp_max, temp_mean):
    """Hargreaves-Samani model for ET0 estimation from temperature input.
    
    Params
    ------
    - lat_deg: float
    - dates: timestamp
    - temp_*: float
    
    
    """
    lat = deg2rad(lat_deg)  # Convert latitude in degrees to radians
    day_of_year = dates.dayofyear
    sol_decli = sol_dec(day_of_year) # Solar declination
    sha = sunset_hour_angle(lat, sol_decli)
    ird = inv_rel_dist_earth_sun(day_of_year)
    et_radia = et_rad(lat, sol_decli, sha, ird) # Extraterrestrial radiation
    return hargreaves(temp_min, temp_max, temp_mean, et_radia)


## Example use with pandas.DataFrame df
## Final dataframe is made by custom function timeseries

# lat_deg = 44.570842547510622 # latitude of Budrio (deg)
# temp_min = df.min()['Temperatura[°C]'].values
# temp_max = df.max()['Temperatura[°C]'].values
# temp_mean = df.mean()['Temperatura[°C]'].values
# dates = df.asfreq().index
# eto = timeseries( dates,
#                  [ hargre(lat_deg, dates[i] , temp_min[i], temp_max[i], temp_mean[i])
#                   for i in range(len(dates)) ] )
# eto_df = pd.DataFrame(eto).rename(columns={0:'Date',1:'EPOT'}).set_index('Date')