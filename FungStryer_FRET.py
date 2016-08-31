# Module containing functions for Fung and Stryer numerical calculations of 2D FRET

# Written by Eric Senning (2016). Last modified August 2016.
# FungStryer_FRET.py is licensed under a 
# Creative Commons Attribution-ShareAlike 3.0 Unported License.

# Import statement
import numpy as np
import scipy as sp
from scipy import integrate

# Fluorescence intensity decay of donor in absence of acceptor:
def FLUORd(t,tau):
	F=np.exp(-1.0*(t/tau))
	return F
	
# Fluorescence intensity decay of donor in presence of acceptor:
def FLUORda(Ro,L,C,t,tau):
	ettSint= lambda x: 2*sp.pi*x*(1.0-np.exp(-(t/tau)*(Ro/x)**6))
	S,err=integrate.quad(ettSint,L,np.inf)
	F=np.exp(-t/tau)*np.exp(-C*S)
	return F

# Calculation of FRET in 2 dimensions for a given distance of closest approach (L)	
def FS_FRET2D_vs_AccDens(Ro,L,accDensities,kvalues,kint,tautime):
	n_accD=np.size(accDensities)
	Ef=[]
	for accD in accDensities:
		sumFda=0.0
		sumFd=0.0
		for kval in kvalues:
			sumFda=kint*FLUORda(Ro,L,accD,kval,tautime)+sumFda
			sumFd=kint*FLUORd(kval,tautime)+sumFd
		energyTransfer=1.0-(sumFda/sumFd)
		Ef.append(energyTransfer)
	return Ef
	
