from __future__ import division
import matplotlib.pyplot as plt
from auxiliary import dspace
import numpy as np
import pickle

class TTsolution(object):
    def __init__(self, scan, refl, scantype, crystal, hkl, R_bend, thickness, E0_or_th):
        self._scan = scan #tuple
        self._reflectivity = refl
        if scantype == 'energy' or scantype == 'wavelength' or scantype == 'angle': 
            self.scantype = scantype
        else:
            raise TypeError()
        self.crystal = crystal.lower()
        self.hkl = hkl
        self._R_bend = R_bend #tuple       
        self._thickness = thickness #tuple

        hc=1.23984193*0.001
        d=dspace(hkl,crystal)*1e-10

        if scantype == 'energy':
            self._th = E0_or_th
            E0=hc/(2*d*np.sin(np.radians(E0_or_th[0])))*1e-6
            wavelength = 2*d*np.sin(np.radians(E0_or_th[0]))

            self._E0 = (E0, 'keV')
            self._lambda0 = (wavelength, 'nm')
        elif scantype == 'wavelength':
            self._th = E0_or_th
            E0=hc/(2*d*np.sin(np.radians(E0_or_th[0])))*1e-6
            wavelength = 2*d*np.sin(np.radians(E0_or_th[0]))

            self._E0 = (E0, 'keV')
            self._lambda0 = (wavelength*1e9, 'nm')
        else:
            self._E0 = E0_or_th
            self._lambda0 = (hc/self._E0[0]*1e9, 'nm')
            self._th = (np.degrees(np.arcsin(hc/(2*d*E0_or_th[0])*1e-6)),'deg')            
   
    @property
    def scan(self):
        return self._scan[0]

    @property
    def reflectivity(self):
        return self._reflectivity

    @property
    def R_bend(self):
        return self._R_bend[0]

    @property
    def thickness(self):
        return self._thickness[0]

    @property
    def E0(self):
        return self._E0[0]

    @property
    def lambda0(self):
        return self._lambda0[0]

    @property
    def th(self):
        return self._th[0]

    def plot(self):
        plt.plot(self._scan[0],self._reflectivity)
        plt.xlabel(self.scantype + ' (' + self._scan[1] + ')')
        plt.ylabel('reflectivity')
        if self.scantype == 'energy' or self.scantype == 'wavelength':
            at_string = str(self._th[0]) + ' ' + self._th[1]     
        else:
            at_string = str(self._E0[0]) + ' ' + self._E0[1] 

        plt.title('Reflectivity of '+ self.crystal + str(self.hkl) + ' @ '+at_string)
        plt.show(block = False)

def save_solution(solution,filename):
    savefile = open(filename,'w')
    pickle.dump(solution,savefile)
    savefile.close()

def load_solution(filename):
    savefile = open(filename,'r')
    solution = pickle.load(savefile)
    savefile.close()
    return solution
