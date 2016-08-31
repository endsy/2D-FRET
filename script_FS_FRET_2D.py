# Script for Fung and Stryer 2D FRET (Fung and Stryer, Biochemistry (1978))
# with acceptor density as abscissa. Parameters accepted in the numerical calculation
# of FRET are, Forster radius (Ro), list of acceptor densities (accDensities), distance of 
# closest approach (L), donor fluorescence lifetime (tautime), time range of integration (kvalues),
# interval of time integration (kint).

# Written by Eric Senning (2016). Last modified August 2016.
# script_FS_FRET_2D.py is licensed under a 
# Creative Commons Attribution-ShareAlike 3.0 Unported License.

#import statement
import scipy as sp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy.matlib import rand,zeros,ones,empty,eye
# import sys
# print sys.path
import FungStryer_FRET as FRET

# Set variables for numerical calculation of FRET
Ro=80.0 #Forster Radius
accDensities=np.linspace(0,.00016,32) #list of acceptor densities in Ang^-2
L=40.0 #distance of closest approach in Ang
tautime=2.0 # this is the lifetime of the donor in nanoseconds
kvalues=np.linspace(0,10,101)# range of integration in nanoseconds
kint=kvalues[1]-kvalues[0]# interval in integration range
# print accDensities
# print "The interval of integration is: %r" % (kint)
# print "The closest approach is set to %r" % (L)
print np.size(accDensities)

# Calculate FRET
EneTra40=FRET.FS_FRET2D_vs_AccDens(Ro,L,accDensities,kvalues,kint,tautime)
L=60
EneTra60=FRET.FS_FRET2D_vs_AccDens(Ro,L,accDensities,kvalues,kint,tautime)
L=100
EneTra100=FRET.FS_FRET2D_vs_AccDens(Ro,L,accDensities,kvalues,kint,tautime)

# print accDensities
# print "---------------"
# print EneTra

#Plot the FRET versus acceptor density curves.
plt.figure(1)
accDum2=accDensities*(10000)**2
plt.plot(accDum2,EneTra40,'r-',accDum2,EneTra60,'m-',accDum2,EneTra100,'b-')
plt.xlabel('Acceptor Density (um-2)')
plt.ylabel('Energy Transfer')
plt.title('Energy transfer versus acceptor density')
plt.text(1000, .4, r'Ro = 50 L =100')
plt.text(1000, .6, r'Ro = 50 L =60')
plt.text(1000, .8, r'Ro = 50 L =40')
plt.axis([0, 18000, 0, 0.9])
plt.show()


