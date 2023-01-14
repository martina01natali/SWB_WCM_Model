# EPOT
def hamon(tavg, jdate, lat, par=1.2):
    # ' @title Hamon Potential Evapotranspiration Equation
    # ' @description The Hamon method is also considered as one of the simplest estimates
    # '     of potential Evapotranspiration.
    # ' @param par proportionality coefficient (unitless),
    # PAR=1.2 as in https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1752-1688.2005.tb03759.x
    # ' @param tavg vector of mean daily temperature (deg C)
    # ' @param lat latitude ()
    # ' @param jdate a day number of the year (julian day of the year)  timetj = df.index.day
    # ' @return outputs potential evapotranspiration (mm day-1)
    # ' @details For details see Haith and Shoemaker (1987) 

    var_theta = 0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.0086 * (jdate - 186)))
    var_pi = np.arcsin(0.39795 * np.cos(var_theta))
    daylighthr = 24 - 24 / np.pi * np.arccos((np.sin(0.8333 * np.pi / 180) + np.sin(lat * np.pi / 180) *
                                              np.sin(var_pi)) / (np.cos(lat * np.pi / 180) * np.cos(var_pi)))
    esat = 0.611 * np.exp(17.27 * tavg / (237.3 + tavg))
    return par * 29.8 * daylighthr * (esat / (tavg + 273.2))