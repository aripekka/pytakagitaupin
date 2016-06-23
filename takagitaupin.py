from __future__ import division
import numpy as np
import os
from auxiliary import *
from scipy.integrate import odeint, ode
import matplotlib.pyplot as plt
from smatrix import compute_S_matrix, compute_S_matrix_fast, mean_poisson
from ttsolution import TTsolution
import sys

def takagitaupin(scantype,scan,constant,hkl,crystal,thickness,bending = 'None'):
    '''
    1D TT-solver.
    
    Input:
    scantype = 'energy' or 'angle'
    scan =  relative to the Bragg's energy in meV (energy scan) OR relative to the Bragg's angle in arcsec (angle scan)
    constant = incidence angle in degrees (energy scan) OR photon energy in keV (angle scan) 
    hkl = [h,k,l] (Miller indices)
    crystal = currently only 'si' is supported
    thickness = crystal thickness in microns
    bending = 'None' OR ('spherical',R_bend) OR ('cylindrical',R1,R2), where R_bend, R1, and R2 are in meters 
    '''

    if scantype == 'energy':
        is_escan = True
        scantype = 'energy'
    elif scantype == 'angle' or scantype == 'angular':
        is_escan = False
        scantype = 'angle'

    #type conversions
    scan=np.array(scan)
        
    #Unit conversions
    thickness_tuple = (thickness, 'microns')
    thickness = thickness*1e-6 #wafer thickness in meters

    #constants
    hc=1.23984193*0.001 #in meV/m
    d=dspace(hkl,crystal)*1e-10 #in m

    #Setup scan variables and constants
    if is_escan:
        escan=scan

        th=np.radians(constant)

        #Direction cosines
        gamma0=np.sin(th)
        gammah=-np.sin(th)

        #Conversion of incident photon energy to wavelength
        E0 = hc/(2*d*np.sin(th)) #in meV

        wavelength = hc/(E0+escan) #in m
    else:
        E0 = constant*1e6 #in meV

        wavelength = hc/E0 #in m

        if not hc/(2*d*E0) > 1:
            th = np.arcsin(hc/(2*d*E0))
        else:
            th = np.pi/2

        ascan = scan*np.pi/648000 #from arcsec to rad

        #Direction cosines
        gamma0=np.sin(th+ascan)
        gammah=-np.sin(th+ascan)

    #construct the path for chitables    
    hklstring = str(hkl[0]) + '_' + str(hkl[1]) + '_' + str(hkl[2])
    filename = 'chitable_' + crystal.lower() + '_' + hklstring + '.dat'
    filestring = os.path.join(os.path.dirname(__file__),'chitables_300K',filename)

    #load the chitable
    try:
        chi = np.loadtxt(filestring)
    except:
        print 'Error loading chitable! Check that ' + filestring \
              + ' exists and is correctly formatted!'
        raise Exception()        

    #conversion to meV
    chienergy = chi[:,0]*1e6

    print 'Computing elastic line for ' + str(hkl) + ' reflection of ' \
          + crystal[0].upper() + crystal[1:].lower() + '-crystal'
    if is_escan:
        print 'Scanning the incident energy'  
    else:
        print 'Scanning the incidence angle'  

    #Polarization (TODO: include pi-polarization)
    C = 1;
    print 'Assuming sigma-polarization'

    #Interpolation
    if is_escan:
        chi0 = np.interp(E0+escan, chienergy, chi[:,1]) + 1j*np.interp(E0+escan, chienergy, chi[:,2])
        chih = np.interp(E0+escan, chienergy, chi[:,3]) + 1j*np.interp(E0+escan, chienergy, chi[:,4])
        chihbar = np.interp(E0+escan, chienergy, chi[:,5]) + 1j*np.interp(E0+escan, chienergy, chi[:,6])
    else:
        chi0 = np.interp(E0, chienergy, chi[:,1]) + 1j*np.interp(E0, chienergy, chi[:,2])
        chih = np.interp(E0, chienergy, chi[:,3]) + 1j*np.interp(E0, chienergy, chi[:,4])
        chihbar = np.interp(E0, chienergy, chi[:,5]) + 1j*np.interp(E0, chienergy, chi[:,6])

    #Deviation from backscattering
    deltawavelength = wavelength-2*d
    if is_escan:
        th2 = th
    else:
        th2 = th+ascan
    #if is_escan:
    #    deltath = th-np.pi/2
    #else:
    #    deltath = th+ascan-np.pi/2

    #Extinction length
    L = wavelength * np.sqrt(gamma0*np.abs(gammah)) / (np.abs(C)*np.sqrt(chih*chihbar))

    #Incidence parameter
    eta = np.sqrt(gamma0/np.abs(gammah)) / (np.abs(C)*np.sqrt(chih*chihbar)) \
        * (-wavelength/d*(wavelength/(2*d)-np.sin(th2)) - chi0*(gammah/gamma0-1)/2)
    #eta = np.sqrt(gamma0/np.abs(gammah)) / (np.abs(C)*np.sqrt(chih*chihbar)) \
    #    * (-wavelength/d*(deltawavelength/(2*d)+deltath**2/2-deltath**4/24+deltath**6/720) - chi0*(gammah/gamma0-1)/2)

    #normalization coefficient
    normcoef = np.sqrt(chih*chihbar)/chihbar*np.sign(C)*np.sqrt(gamma0/np.abs(gammah))

    #Calculate mean poisson's ratio
    nu = 0

    if not bending == 'None':
        #TODO: different bendings have their own quirks, check for cylindrical
        S_matrix, C_matrix = compute_S_matrix_fast(hkl,crystal)
        #nu = mean_poisson(S_matrix)
        
        #test
        S=S_matrix
        
        if bending[0] == 'cylindrical':
            if bending[1] == 'inf':
                invR1 = 0
            else:
                invR1 = 1/bending[1]

            if bending[2] == 'inf':
                invR2 = 0
            else:            
                invR2 = 1/bending[2]            
        
        elif bending[0] == 'spherical':
            if bending[1] == 'inf':
                invR1 = 0
                invR2 = 0
            else:
                invR1 = 1/bending[1]
                invR2 = 1/bending[1]
        
        #Parameter according to http://arxiv.org/abs/1502.03059
        bending_parameter = S[2,0]*(S[0,1]*invR2-S[1,1]*invR1)+S[2,1]*(S[1,0]*invR1-S[0,0]*invR2)
        bending_parameter = -0.5*bending_parameter/(S[0,1]*S[1,0]-S[0,0]*S[1,1])
        print bending_parameter

    #INTEGRATION
    reflectivity=[]

    #Define ODE and its Jacobian
    def tt_equation(z,ksi,L,gamma0,gammah,eta,d,bending,thickness,nu):
        if bending == 'None':
            return np.pi*1j/L*(ksi**2-2*(np.sign(gammah)*eta)*ksi-np.sign(gammah))
        else:
            return np.pi*1j/L*(ksi**2-2*(np.sign(gammah)*eta+L*2*bending_parameter*(z-thickness/2)/d)*ksi-np.sign(gammah))


    def tt_jacobian(z,ksi,L,gamma0,gammah,eta,d,bending,thickness,nu):
        if bending == 'None':
            return np.pi*1j/L*(2*ksi-2*(np.sign(gammah)*eta))
        else:
            return np.pi*1j/L*(2*ksi-2*(np.sign(gammah)*eta+L*2*bending_parameter*(z-thickness/2)/d))


    #Solve the equation

    sys.stdout.write('Solving...0%')
    sys.stdout.flush()
    
    for step in xrange(len(scan)):
        def tt2solve(z,ksi):
            if is_escan:
                return tt_equation(z,ksi,L[step],gamma0,gammah,eta[step],d,bending,thickness,nu)
            else:
                return tt_equation(z,ksi,L[step],gamma0[step],gammah[step],eta[step],d,bending,thickness,nu)

        def jac(z,ksi):
            if is_escan:
                return tt_jacobian(z,ksi,L[step],gamma0,gammah,eta[step],d,bending,thickness,nu)
            else:
                return tt_jacobian(z,ksi,L[step],gamma0[step],gammah[step],eta[step],d,bending,thickness,nu)

        r=ode(tt2solve,jac).set_integrator('zvode',method='bdf',with_jacobian=True,min_step=1e-10,max_step=1e-4,nsteps=50000)
        r.set_initial_value(0,thickness)
        res=r.integrate(0)     
   
        reflectivity.append(np.abs(normcoef[step]*res[0])**2)

        sys.stdout.write('\rSolving...%0.1f%%' % (100*(step+1)/len(scan),))  
        sys.stdout.flush()

    sys.stdout.write('\r\nDone.\n')
    sys.stdout.flush()

    #solution class    
    if is_escan:
        scan = (scan, 'meV')
        constant = (constant,'deg')
    else:
        scan = (scan, 'arcsec')
        constant = (constant,'keV')

    #TODO: add also the type of bending to the ttsolution
    if bending == 'None':
        R_bend = 0
    else:
        R_bend = bending[1]
    
    result = TTsolution(scan,reflectivity,scantype,crystal.lower(),hkl,(R_bend,'m'),thickness_tuple,constant)

    return result
