# -*- coding: utf-8 -*-
"""
Created on Tue May 18 08:27:27 2021

@author: Aaron
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.interpolate import CubicSpline
# from scipy.interpolate import interp1d

from color_chart_funcs import refspec, refrgb

#interpolation functions
from color_chart_data import csn, cso, cssn, cssk, air, pfluor, csciex, csciey, csciez

lam_step = 1
lam = np.arange(380,781,lam_step)

nnit = csn(lam)
nox = cso(lam)
nsi = cssn(lam)
ksi = cssk(lam)
nair = air(lam)

p0 = pfluor(lam)

ciex = csciex(lam)
ciey = csciey(lam)
ciez = csciez(lam)

dt=1
nitride_thicknesses = np.arange(0,221,dt)
rgb_array1 = np.zeros((nitride_thicknesses.size,3))
rgb_array2 = np.zeros((nitride_thicknesses.size,3))

for i, d in enumerate(nitride_thicknesses):
    # rgb_array1[i,:] = refrgb(lam,p0,nair,0,'mixed',nsi-1j*ksi,[nnit,d])
    rgb_array1[i,:] = refrgb(lam,p0,nair,0,None,nsi-1j*ksi,[nnit,d])
    rgb_array2[i,:] = refrgb(lam,p0,nair,0,None,nair,[nnit,d])


fig, ax = plt.subplots(figsize=(10,4))

dif = lam_step
yf1 = np.asarray([0, 0, 200, 200])
yf2 = np.asarray([200, 200, 400, 400])
for i, d in enumerate(nitride_thicknesses):
    xf = np.array([d-dt/2,d+dt/2,d+dt/2,d-dt/2])
    ax.fill(xf,yf1,facecolor=rgb_array1[i,:],edgecolor='none')
    ax.fill(xf,yf2,facecolor=rgb_array2[i,:],edgecolor='none')
    
ax.set_xticks(np.arange(0,221,20))
ax.set_xlabel('Nitride Thickness (nm)')
ax.set_yticks([100,300])
ax.set_yticklabels(['Air/SiN/Si','Air/SiN/Air'])
ax.grid(axis='x', color='w', linestyle=(0, (5.0, 40.0)), linewidth=0.5)
    
# plt.savefig('sin_and_sin-on-si.png', dpi=600)