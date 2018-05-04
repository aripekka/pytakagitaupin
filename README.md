# pytakagitaupin

##THIS PROJECT IS OBSOLETE. FOR NEW VERSION, SEE https://github.com/aripekka/pyTTE

Python solver for 1D Takagi-Taupin equation

Usage:

The function 'takagitaupin' in takagitaupin.py solves the reflectivity of a crystal as a function of the photon energy or the incidence angle. The crystal can be undeformed or spherically bent. Currently only symmetric Bragg case with sigma-polarisation is supported

To make the energy scan, the takagitaupin is used as follows:

sol = takagitaupin('energy',energy-vector,incidence-angle,[h,k,l],crystal-string,bending-radius,crystal-thickness)

energy-vector (in meV)
incidence-angle (in degrees)
[h,k,l] = the Miller indices (note that the reflection might have to be added to chitables)
crystal string = currently 'si' or 'ge'
bending-radius (spherical bendin in meters or 'inf' for undeformed)
crystal-thickness (in microns)

For an example, to compute the reflectivity of a spherically bent Si(660) analyser (with Rb = 1 m) using only the 
depth-dependent strain near back-scattering conditions (thetab = 88.4 deg), the function call is

sol = takagitaupin('energy',linspace(-800,200,200),88.4,[6,6,0],'si',1,300)

Similarly for the angular scan:

sol = takagitaupin('angle',angle-vector,photon_energy,[h,k,l],crystal-string,bending-radius,crystal-thickness)

Here angle-vector is in arcsecs and photon_energy in keV.

sol is a TTsolution object that can be inspected quickly with sol.plot(). 
The scan and reflectivity vectors are obtained respectively via properties sol.scan and sol.reflectivity.



