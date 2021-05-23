# -*- coding: utf-8 -*-
"""
Created on Mon May 17 23:09:08 2021

@author: Aaron
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d


file_nitride = 'nitride-Philipp.csv' # wl n
file_oxide = 'oxide-Gao.csv' # wl n
file_silicon = 'silicon-Vuye-20C.csv' # wl n k

file_fluorescent = 'standard-fluorescent-normalized-digitized.csv'
file_ciexyz = 'ciexyzj-mod.csv'
file_tristimulus = 'CIE1931_tristimulus_5nm_edit.csv'


nitride = np.genfromtxt(file_nitride, dtype=float, delimiter=',', skip_header=1)
csn = CubicSpline(nitride[:,0]*1000,nitride[:,1]) # for wavelength in nm

oxide = np.genfromtxt(file_oxide, dtype=float, delimiter=',', skip_header=1)
cso = CubicSpline(oxide[:,0]*1000, oxide[:,1]) # for wavelength in nm

silicon = np.genfromtxt(file_silicon, dtype=float, delimiter=',', skip_header=1)
cssn = CubicSpline(silicon[:,0]*1000,silicon[:,1]) # n for wavelength in nm
cssk = CubicSpline(silicon[:,0]*1000,silicon[:,2]) # k for wavelength in nm

fluor = np.genfromtxt(file_fluorescent, dtype=float, delimiter=',')
pfluor = interp1d(fluor[:,0],fluor[:,1])
# csfluor = CubicSpline(fluor[:,0], fluor[:,1]) # normalized incident power

# ciexyz = np.genfromtxt(file_ciexyz, dtype=float, delimiter=',')
ciexyz = np.genfromtxt(file_tristimulus, dtype=float, delimiter = ',',skip_header=1)
csciex = CubicSpline(ciexyz[:,0],ciexyz[:,1])
csciey = CubicSpline(ciexyz[:,0],ciexyz[:,2])
csciez = CubicSpline(ciexyz[:,0],ciexyz[:,3])

matrix_srgb_fluor = np.asarray([[23.42, -9.83, 1.16],[-11.11, 19.02, -4.24],[-3.6, 0.42, 21.96]])
# xyz_to_rgb = np.linalg.inv(matrix_srgb_fluor) # perhaps there's a mistake in the paper
xyz_to_rgb = matrix_srgb_fluor

lam_step = 1

lam = np.arange(380,781,lam_step)
nnit = csn(lam)
nox = cso(lam)
nsi = cssn(lam)
ksi = cssk(lam)
nair = np.ones(lam.shape)
def air(lam):
    return np.ones(lam.shape)

p0 = pfluor(lam)

ciex = csciex(lam)
ciey = csciey(lam)
ciez = csciez(lam)


