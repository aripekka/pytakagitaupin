from __future__ import division
import numpy as np
from auxiliary import rotation_matrix2 as rotation_matrix

def compute_S_matrix_fast(zdir, xtal):
    '''
        Computes the compliance and stiffness matrices S and C a given z-direction.
        The x- and y-directions are determined automatically 
    '''
    xtal = xtal.lower()

    #TODO: read the elastic constants from a file 
    if xtal == 'ge':
    	c1111, c1122, c2323 = 1.292, 0.479, 0.670
    elif xtal == 'si':
    	c1111, c1122, c2323 = 1.657, 0.639, 0.796
    else:
        raise ValueError('Elastic parameters for the crystal not found!')

    #TODO: generalize to other systems alongside the cubic as well
    Cc = np.zeros((3,3,3,3))

    Cc[0,0,0,0], Cc[1,1,1,1], Cc[2,2,2,2] = c1111, c1111, c1111
    Cc[0,0,1,1], Cc[0,0,2,2], Cc[1,1,0,0] = c1122, c1122, c1122
    Cc[1,1,2,2], Cc[2,2,0,0], Cc[2,2,1,1] = c1122, c1122, c1122

    Cc[0,2,0,2], Cc[2,0,0,2], Cc[0,2,2,0], Cc[2,0,2,0] = c2323, c2323, c2323, c2323
    Cc[1,2,1,2], Cc[2,1,1,2], Cc[1,2,2,1], Cc[2,1,2,1] = c2323, c2323, c2323, c2323
    Cc[0,1,0,1], Cc[1,0,0,1], Cc[0,1,1,0], Cc[1,0,1,0] = c2323, c2323, c2323, c2323

    Q = rotation_matrix(zdir)
   
    #Rotate the tensor

    Crot = np.zeros((3,3,3,3))
    for i in range(3):
     for j in range(3):
      for k in range(3):
       for l in range(3):
        for p in range(3):
         for q in range(3):
          for r in range(3):
           for s in range(3):
            Crot[i,j,k,l]=Crot[i,j,k,l]+Q[i,p]*Q[j,q]*Q[k,r]*Q[l,s]*Cc[p,q,r,s]     

    C = np.array([
        [Crot[0,0,0,0], Crot[0,0,1,1], Crot[0,0,2,2], Crot[0,0,1,2], Crot[0,0,0,2], Crot[0,0,0,1]],
        [Crot[1,1,0,0], Crot[1,1,1,1], Crot[1,1,2,2], Crot[1,1,1,2], Crot[1,1,0,2], Crot[1,1,0,1]],
        [Crot[2,2,0,0], Crot[2,2,1,1], Crot[2,2,2,2], Crot[2,2,1,2], Crot[2,2,0,2], Crot[2,2,0,1]],
        [Crot[2,1,0,0], Crot[2,1,1,1], Crot[2,1,2,2], Crot[1,2,1,2], Crot[1,2,0,2], Crot[1,2,0,1]],
        [Crot[2,0,0,0], Crot[2,0,1,1], Crot[2,0,2,2], Crot[2,0,1,2], Crot[0,2,0,2], Crot[2,0,0,1]],
        [Crot[1,0,0,0], Crot[1,0,1,1], Crot[1,0,2,2], Crot[1,0,1,2], Crot[1,0,0,2], Crot[0,1,0,1]]
    ])

    C=C*1e11 #in pascal

    S = np.linalg.inv(C)

    return S, C

def mean_poisson(S):
    '''Computes the mean of the Poisson ratio over 2*pi (xy-plane) from given S matrix'''

    #Rotated matrix components (computed with Maxima)
    x = np.linspace(0,2*np.pi,500)

    Sr11 = S[0,0]*np.cos(x)**4 + S[1,1]*np.sin(x)**4 \
         + 2*np.sin(x)*np.cos(x)*(S[0,5]*np.cos(x)**2 + S[1,5]*np.sin(x)**2) \
         + (2*S[0,1]+S[5,5])*np.cos(x)**2*np.sin(x)**2
    Sr31 = S[2,0]*np.cos(x)**2 + S[2,1]*np.sin(x)**2 + S[2,5]*np.cos(x)*np.sin(x)


    return np.trapz(-Sr31/Sr11,x)/(2*np.pi)

def compute_S_matrix(zdir,xtal):
    if xtal.lower() == 'ge':
    	c1111, c1122, c2323 = 1.292, 0.479, 0.670
    elif xtal.lower() == 'si':
    	c1111, c1122, c2323 = 1.657, 0.639, 0.796        

    Cc = np.zeros((3,3,3,3))

    Cc[0,0,0,0], Cc[1,1,1,1], Cc[2,2,2,2] = c1111, c1111, c1111
    Cc[0,0,1,1], Cc[0,0,2,2], Cc[1,1,0,0] = c1122, c1122, c1122
    Cc[1,1,2,2], Cc[2,2,0,0], Cc[2,2,1,1] = c1122, c1122, c1122

    Cc[0,2,0,2], Cc[2,0,0,2], Cc[0,2,2,0], Cc[2,0,2,0] = c2323, c2323, c2323, c2323
    Cc[1,2,1,2], Cc[2,1,1,2], Cc[1,2,2,1], Cc[2,1,2,1] = c2323, c2323, c2323, c2323
    Cc[0,1,0,1], Cc[1,0,0,1], Cc[0,1,1,0], Cc[1,0,1,0] = c2323, c2323, c2323, c2323

    N = 500
    theta = np.linspace(0,2*np.pi,N)
    Q2 = rotation_matrix(zdir)

    S,C=[],[]

    for th in theta:
        Q = np.array([[np.cos(th), np.sin(th), 0],[-np.sin(th), np.cos(th), 0], [0, 0, 1]])
        Q = np.dot(Q, Q2)

        Crot = np.zeros((3,3,3,3))
        for i in range(3):
         for j in range(3):
          for k in range(3):
           for l in range(3):
            for p in range(3):
             for q in range(3):
              for r in range(3):
               for s in range(3):
                Crot[i,j,k,l]=Crot[i,j,k,l]+Q[i,p]*Q[j,q]*Q[k,r]*Q[l,s]*Cc[p,q,r,s]     

        Cr = np.array([
            [Crot[0,0,0,0], Crot[0,0,1,1], Crot[0,0,2,2], Crot[0,0,1,2], Crot[0,0,0,2], Crot[0,0,0,1]],
            [Crot[1,1,0,0], Crot[1,1,1,1], Crot[1,1,2,2], Crot[1,1,1,2], Crot[1,1,0,2], Crot[1,1,0,1]],
            [Crot[2,2,0,0], Crot[2,2,1,1], Crot[2,2,2,2], Crot[2,2,1,2], Crot[2,2,0,2], Crot[2,2,0,1]],
            [Crot[2,1,0,0], Crot[2,1,1,1], Crot[2,1,2,2], Crot[1,2,1,2], Crot[1,2,0,2], Crot[1,2,0,1]],
            [Crot[2,0,0,0], Crot[2,0,1,1], Crot[2,0,2,2], Crot[2,0,1,2], Crot[0,2,0,2], Crot[2,0,0,1]],
            [Crot[1,0,0,0], Crot[1,0,1,1], Crot[1,0,2,2], Crot[1,0,1,2], Crot[1,0,0,2], Crot[0,1,0,1]]
        ])
    
        Cr=Cr*1e11 #in pascal

        #By definitions
        S.append(np.linalg.inv(Cr))
        C.append(Cr)

    angle = 180*theta/np.pi

    return angle, S, C
