from __future__ import division
from smatrix import compute_S_matrix
import numpy as np
import matplotlib.pyplot as plt

def analyser_elastic(ttsolution,diameter):
    '''
    Code to compute the elastic line of a spherically bent analyser (diameter in mm). 
    '''
    #TODO generalise to angular scans
    if not ttsolution.scantype == 'energy':
        raise ValueError('Currently the code supports only energy scans')

    #Takagi-Taupin curve
    scan = ttsolution.scan*1e-3 #scan in eV
    refl = ttsolution.reflectivity

    #Create sampling grid
    radius = 0.5*diameter*1e-3
    x = np.linspace(-radius,radius,250)
    X, Y = np.meshgrid(x,x)

    #conversion to polar coordinates
    R = np.sqrt(X**2 + Y**2)
    TH = np.arctan2(Y, X)    

    S=compute_S_matrix_fast(ttsolution.hkl,ttsolution.crystal)[0]

    #the strain field
    Rb = ttsolution.R_bend
    denominator = 5*(S[0,0]+S[1,1])+6*S[1,0]+S[5,5]
    A = -4/3*(S[2,0]+S[2,1])/denominator
    B = -2/3*np.sqrt((S[2,1]-S[2,0])**2 + S[2,5]**2)/denominator
    phi0 = np.arctan(S[2,5]/(S[2,1]-S[2,0]))

    epsilon = (A+B*np.cos(2*TH+phi0))*(R/Rb)**2

    #energy shift map
    deltaE = -epsilon*ttsolution.E0*1e3 

    #Collect energy shifts belonging to the active surface of the analyser
    indices = np.nonzero(R < radius)
    deltaE = deltaE[indices]

    #grid for for the elastic line
    energy = np.linspace(deltaE.min()+scan.min(),deltaE.max()+scan.max(),250)
    I = np.zeros(energy.size)

    for dE in deltaE:
        I = I + np.interp(energy-dE,scan,refl,left=0,right=0)
        
    #normalization
    I = I/len(deltaE)*(diameter/100)**2

    #plt.plot(energy,I)
    #plt.show(block = False)

    return energy,I
