"""
c.massari - s.modanesi 08/11/2022
define WCM(PAR[i], data_in) 
#può essere usato con modello di suolo + WCM o WCM o modello di suolo
con
"""
import pyswarms as ps
import hydroeval as he
import numpy as np
import matplotlib.pyplot as plt



def WCM(PAR, data_in):

    A,B,C,D=PAR
    SM, LAI,t_deg,obs=data_in
    theta = t_deg * np.pi / 180.
    sig0s_dB = C + D * SM #define bare soil backscatter
    sig0s=10**(sig0s_dB/10) #from dB to linear scale
    T2 = np.exp((-2 * B * LAI) / np.cos(theta)) #two-way attenuation from the vegetation layer
    sig0v = A * LAI * np.cos(theta) * (1 - T2) #define backscatter from the vegetation
    sig0_lin = T2 * sig0s + sig0v
    sig0=10*np.log10(sig0_lin) #from linear scale to dB
    
    OUT=he.evaluator(he.kge, sig0, obs) #OUT is kge, r, alpha, betha
    KGE=OUT[0,:];

    return KGE

#----------------TRY an EXAMPLE with random timeseries
PAR=[0.4,0.4,-20,40] #some guess values of WCM Parameters
SM=np.random.uniform(low=0.05, high=0.5, size=(100,))
LAI=np.random.uniform(low=0, high=6, size=(100,))
t_deg = np.random.uniform(low=37, high=37, size=(100,))
obs = np.random.uniform(low=-20, high=-5, size=(100,))
plt.plot(SM)
plt.plot(obs)
plt.plot(t_deg)

data_in=[SM, LAI,t_deg,obs]
KGE = WCM(PAR, data_in)

#--------------START calibration (here using random timeseries)
bnds1 = (np.array([0.4, 0.4, -10,80]),
                 np.array([0, 0, -35, 15]))

options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}

def func0(PAR):
    global data_in
    n_particles = PAR.shape[0]
    err = np.zeros(n_particles)
    for i in range(n_particles):
        KGE = WCM(PAR[i], data_in)        
        err[i] = 1 - KGE
    return err


optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=4, options=options, bounds=bnds1)

cost, PARn = optimizer.optimize(func0, 10) #revise

#------------- model run validation
results = WCM(PARn, data_in)

