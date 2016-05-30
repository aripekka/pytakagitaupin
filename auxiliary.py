from __future__ import division
import numpy as np

def dspace(hkl,xtal):
    constants={'si' : 5.43102088}
    divisor=np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2)
    return constants[xtal.lower()]/divisor

def rotation_matrix(hkl):
    ##Calculates the rotation matrix R that aligns the desired crystal direction parallel to z-axis.
    ##NOTE: This version works only for the cubic systems!

    if hkl[0]:
	    yx = hkl[1]/hkl[0]
	    R_z = np.array([[1, yx, 0], [-yx, 1, 0], [0, 0, np.sqrt(1+yx**2)]])/np.sqrt(1+yx**2)
    elif hkl[1]>0:
	    R_z = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    elif hkl[1]<0:
	    R_z = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]);
    else:
	    R_z = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])

    zr=hkl[2]/np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2)

    R_y=np.array([[zr, 0, -np.sqrt(1-zr**2)], [0, 1, 0], [np.sqrt(1-zr**2), 0, zr]]);

    return np.dot(R_y,R_z)

def rotation_matrix2(hkl):
    if hkl[0] or hkl[1]:
        #rotation axis
        u = -np.array([[hkl[1]],[-hkl[0]]])/np.sqrt(hkl[0]**2+hkl[1]**2)
        #rotation angle
        th = np.arccos(hkl[2]/np.sqrt(hkl[0]**2+hkl[1]**2+hkl[2]**2))

        #rotation matrix
        R=np.array([[np.cos(th)+u[0]**2*(1-np.cos(th)), u[0]*u[1]*(1-np.cos(th)), u[1]*np.sin(th)],
           [u[0]*u[1]*(1-np.cos(th)), np.cos(th)+u[1]**2*(1-np.cos(th)), -u[0]*np.sin(th)],
           [-u[1]*np.sin(th), u[0]*np.sin(th), np.cos(th)]])
    else:
        R=np.array([[1,0,0],[0,1,0],[0,0,1]])

    return R.transpose()

def gauss(x,x0,fwhm):
    return np.exp(-np.log(2)*((x-x0)/fwhm*2)**2)

