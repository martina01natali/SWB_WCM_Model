import numpy as np

# Part,Frequency(GHz),a0,a1,a2,b0,b1,b2,c0,c1,c2
coeff = {
    4: 
    {
        'real':[ 2.927,-0.012,-0.001,5.505,0.371,0.062,114.826,-0.389,-0.547],
        'img':[0.004,0.001,0.002,0.951,0.005,-0.01,16.759,0.192,0.29]
    },
    6: 
    {
        'real':[1.993,0.002,0.015,38.086,-0.176,-0.633,10.72,1.256,1.522],
        'img':[-0.123,0.002,0.003,7.502,-0.058,-0.116,2.942,0.452,0.543]
    }
}


def hallikainen(freq:int, sand:float, clay:float, water:np.array):
    """
    coeff: dict('real':[],'img':[])
    """
    global coeff
    
    coeff_r = np.array(coeff[freq]['real'])
    real = coeff_r[0]+coeff_r[1]*sand+coeff_r[2]*clay+\
    (coeff_r[3]+coeff_r[4]*sand+coeff_r[5]*clay)*water+\
    (coeff_r[6]+coeff_r[7]*sand+coeff_r[8]*clay)*(water**2)
        
    coeff_i = np.array(coeff[freq]['img'])
    img = coeff_i[0]+coeff_i[1]*sand+coeff_i[2]*clay+\
    (coeff_i[3]+coeff_i[4]*sand+coeff_i[5]*clay)*water+\
    (coeff_i[6]+coeff_i[7]*sand+coeff_i[8]*clay)*(water**2)
    
    return real, img

def doi(freq:int, sand:float, clay:float, water:np.array, angle:float):
    """
    freq [GHz]
    angle [°]
    
    return depth [m]
    """
    c = 299792458 # m/s
    theta  = angle*np.pi/180. # angle [rad]
    
    real = hallikainen(freq, sand=sand, clay=clay, water=water)[0]
    img = hallikainen(freq, sand=sand, clay=clay, water=water)[1]
    
    depth =  c/(2*np.pi*freq*1e9)*(np.sqrt(real)/img)*np.cos(theta)
    
    return depth