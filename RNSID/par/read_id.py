#!/usr/bin/env python

from __future__ import print_function
import  numpy as  np
import matplotlib as mpl
from matplotlib import pyplot as plt
import h5py 
import sys

filename = sys.argv[1]

with h5py.File(filename,'r') as f:
    K=f.attrs['poly K'][0]
    Gamma = f.attrs['poly Gamma'][0]
    axis_ratio = f.attrs['axis_ratio'][0]
    re = f.attrs['re'][0]
    rhoc = f.attrs['rhoc'][0]
    A_diff = f.attrs['A diff'][0]
    Omega = f.attrs['Omega'][0]
    Omega_e = f.attrs['Omega_e'][0]
    SDIV=f.attrs['SDIV'][0]
    MDIV=f.attrs['MDIV'][0]
    rtype = f.attrs['rotation_type'][0]
    eos_type=f.attrs['eos_type'][0]
    eos_file = f.attrs['eos_file'][0]
    Omega_diff = f['Omega_diff'][:]
    alpha = f['alpha'][:]
    energy = f['energy'][:]
    enthalpy = f['enthalpy'][:]
    gama = f['gama'][:]
    mu = f['mu'][:]
    omega = f['omega'][:]
    pressure = f['pressure'][:]
    rho_potential = f['rho_potential'][:]
    s = f['s_qp'][:]
    v2 = f['velocity_sq'][:]

S,MU = np.meshgrid(s,mu,indexing='ij')
NU = np.sqrt(1 - MU**2)
XC = S*MU
ZC = S*NU

plt.pcolor(XC,ZC,np.log(pressure))
plt.savefig('star_compactified.png',bbox_inches='tight')
plt.clf()

R = -re*S/(S-1)
rmax = 2*re
irmax = np.where(R >= rmax)[0][0]
X = (R*MU)[:irmax]
Y = (R*NU)[:irmax]
LOGP = np.log(pressure)[:irmax-1]

plt.pcolor(X,Y,LOGP)
plt.savefig('star_spherical_polar.png',bbox_inches='tight')
