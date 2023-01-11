# WORK IN PROGRESS
# still defined global variables
# need to implement arguments passing inside the function to make it work regardless fo global parameters
#
# may also want to implement custom optimization algorithm


import os
import re
import time
import math
import numpy as np
import pandas as pd
import datetime as dtt
import hydroeval as he


def pso_calib_irri(PAR):
    """Ausiliary function for PSO optimization"""
    global inputs
    global irri
    n_particles = PAR.shape[0]
    err = np.zeros(n_particles)
    for i in range(n_particles):
        WW,IRR,sig0,KGE = IRR_WCM(PAR[i], inputs, irri)
        err[i] = 1 - KGE
    return err