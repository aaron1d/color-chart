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
from color_chart_data import csn, cso, cssn, cssk, csh, air, pfluor, csciex, csciey, csciez



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

dt = 1
nitride_thicknesses = np.arange(0,200,dt)
rgb_array = np.zeros((nitride_thicknesses.size,3))

for i, d in enumerate(nitride_thicknesses):
    rgb_array[i,:] = refrgb(lam,p0,nair,np.pi/10,'mixed',nsi-1j*ksi,[nnit,d])

plt.close()

fig, ax = plt.subplots(figsize=(10,4))

yf = np.array([0, 0, 400, 400])
for i, d in enumerate(nitride_thicknesses):
    xf = np.array([d-dt/2,d+dt/2,d+dt/2,d-dt/2])
    ax.fill(xf,yf,facecolor=rgb_array[i,:],edgecolor='none')

ax.set_xticks(np.arange(0,201,20))
ax.set_xlabel('Nitride Thickness (nm)')
ax.set_yticks([])
ax.set_yticklabels([])
ax.grid(axis='x', color='w', linestyle='--', linewidth=0.5)
