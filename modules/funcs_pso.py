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

# PSO
from pyswarms.utils.plotters import plot_cost_history
from pyswarms.utils.functions.single_obj import sphere
from IPython.display import Image as Image
from pyswarms.utils.functions.single_obj import sphere
from pyswarms.utils.plotters.formatters import Mesher
from pyswarms.utils.functions import single_obj as fx
from pyswarms.utils.plotters import (plot_cost_history, plot_contour, plot_surface)
from pyswarms.backend.handlers import OptionsHandler
