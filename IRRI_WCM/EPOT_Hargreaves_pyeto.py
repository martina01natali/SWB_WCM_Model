# Ref. https://pyeto.readthedocs.io/en/latest/hargreaves.html

import datetime, pyeto

def hargre(lat_rad, date:datetime.date, temp_min, temp_max, temp_mean)
    lat = pyeto.deg2rad(lat_rad)  # Convert latitude in degrees to radians
    day_of_year = date.timetuple().tm_yday
    sol_dec = pyeto.sol_dec(day_of_year)            # Solar declination
    sha = pyeto.sunset_hour_angle(lat, sol_dec)
    ird = pyeto.inv_rel_dist_earth_sun(day_of_year)
    et_rad = pyeto.et_rad(lat, sol_dec, sha, ird)   # Extraterrestrial radiation
    hargreaves(temp_min, temp_max, temp_mean, et_rad)
