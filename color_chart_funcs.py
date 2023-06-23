# -*- coding: utf-8 -*-
"""
Created on Tue May 18 08:23:55 2021

@author: Aaron
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d

from color_chart_data import ciex, ciey, ciez, xyz_to_rgb

def _delta(lam,n,d,alpha=0):
    '''
    The phase accumulated through a given layer, given λ, n(λ), and α

    Parameters
    ----------
    lam : 1D numpy array of shape (nlam,) 
        wavelength range in nanometers.
    n : 1D numpy array of shape (nlam,)
        refractive index at the values of lam
    d : int or float
        layer thickness in nanometers.
    alpha : int or float
        angle of incidence for a given layer.

    Returns
    -------
    1D numpy array of shape (nlam,)

    '''
    return 2*np.pi/lam*n*d*np.cos(alpha)

def _gamma(n,alpha,pol):
    '''
    The angle and polarization-dependent effect refractive index for a layer

    Parameters
    ----------
    n : float or complex
        refractive index.
    alpha : int or float
        angle of incidence for a given layer.
    pol : string
        'TE' or 'TM', otherwise normal incidence.

    Returns
    -------
    float
        The angle and pol-dependent effect refractive index for a layer.

    '''
    
    if pol=='TE':
        return n*np.cos(alpha)
    elif pol=='TM':
        return n/np.cos(alpha)
    else:
        return n

# For each layer, make an n_lambda x 2 x 2 array

def _mmj(lam,n,d,alpha=0,pol=None):
    '''
    Array of 2x2 optical transfer matrices, length nlam

    Parameters
    ----------
    lam : 1D numpy array of shape (nlam,) 
        wavelength range in nanometers.
    n : 1D numpy array of shape (nlam,)
        refractive index at the values of lam
    d : int or float
        layer thickness in nanometers.
    alpha : int or float
        angle of incidence for a given layer.
    pol : string
        'TE' or 'TM', otherwise normal incidence.

    Returns
    -------
    Numpy array of shape (nlam,2,2)
        Optical transfer matrix for a given layer, at all wavelengths

    '''
    
    # delta(λ,n) and gamma(n) for calculating the matrix 
    dj = _delta(lam,n,d,alpha)
    gj = _gamma(n,alpha,pol)
    
    #array of 2x2 arrays
    return np.asarray([[[np.cos(ddj), 1j*np.sin(ddj)/ggj],
                        [1j*np.sin(ddj)*ggj, np.cos(ddj)]] for ddj,ggj in zip(dj,gj)])


def refspec(lam,n0,alpha0,pol,ns,*ndpairs):
    '''
    Calculates the (Intensity) Reflectance spectrum R(λ) given incident 
    material (n0), incident angle, polarization, substrate material (ns)
    and index-thickness pairs for layers of materials.

    Parameters
    ----------
    lam : 1D numpy array of shape (nlam,) 
        wavelength range in nanometers.
    n0 : 1D numpy array of shape (nlam,)
        refractive index of incident material at the values of lam
    alpha0 : int or float
        angle of incidence for incident layer, default normal incidence
    pol : string
        'TE', 'TM', or 'mixed' (50% each), otherwise normal incidence.
    ns : 1D numpy array of shape (nlam,)
        refractive index of substrate material at the values of lam
    *ndpairs : variable number of 2 element lists of [1d numpy array,int/float ]
       The RI n(λ) and thickness d (in nm) of each layer.

    Returns
    -------
    1d numpy array of shape (nlam,)
        Reflectance spectrum.

    '''
    
    # if mixed polarization, assume 50% TE, 50% TM
    if pol=='mixed':
        refspecTE = refspec(lam,n0,alpha0,'TE',ns,*ndpairs)
        refspecTM = refspec(lam,n0,alpha0,'TM',ns,*ndpairs)
        return (refspecTE + refspecTM)/2
    
    # for nx2x2 array, start with identity matrices
    M = np.asarray([[[1,0],[0,1]] for i in range(lam.size)])
    
    # gamma for incident material
    gg0 = _gamma(n0,alpha0,pol)
    
    # incidence angle for substrate material
    alphas = np.arcsin(n0/ns*np.sin(alpha0))
    
    # gamma for substrate material
    ggs = _gamma(ns,alphas,pol)
    
    # calculate the nx2x2 array for each layer, multiply the 2x2's, for each layer
    for nd in ndpairs:
        nj = nd[0]
        dj = nd[1]
        
        alphaj = np.arcsin(n0/nj*np.sin(alpha0))
        mj = _mmj(lam,nj,dj,alphaj,pol)
        
        M = np.matmul(mj,M)
        
        mlist = list(M)

    # calculate the (amplitude) reflection coeffecient given the total matrix
    ref = [(g0*m[0,0]+gs*g0*m[0,1]-m[1,0]-gs*m[1,1])/(g0*m[0,0]+gs*g0*m[0,1]+m[1,0]+gs*m[1,1])\
               for g0,gs,m in zip(gg0,ggs,mlist)]
    
    # return the intensity reflectance
    return np.abs(ref)**2

# def compref(lam,n0,alpha0,pol,ns,*ndpairs):
#     '''
#     Calculates the complex reflection spectrum r(λ) given incident 
#     material (n0), incident angle, polarization, substrate material (ns)
#     and index-thickness pairs for layers of materials.

#     Parameters
#     ----------
#     lam : 1D numpy array of shape (nlam,) 
#         wavelength range in nanometers.
#     n0 : 1D numpy array of shape (nlam,)
#         refractive index of incident material at the values of lam
#     alpha0 : int or float
#         angle of incidence for incident layer, default normal incidence
#     pol : string
#         'TE', 'TM', or 'mixed' (50% each), otherwise normal incidence.
#     ns : 1D numpy array of shape (nlam,)
#         refractive index of substrate material at the values of lam
#     *ndpairs : variable number of 2 element lists of [1d numpy array,int/float ]
#        The RI n(λ) and thickness d (in nm) of each layer.

#     Returns
#     -------
#     1d numpy array of shape (nlam,)
#         Reflectance spectrum.

#     '''
    
#     # if mixed polarization, assume 50% TE, 50% TM
#     if pol=='mixed':
#         refspecTE = refspec(lam,n0,alpha0,'TE',ns,*ndpairs)
#         refspecTM = refspec(lam,n0,alpha0,'TM',ns,*ndpairs)
#         return (refspecTE + refspecTM)/2
    
#     # for nx2x2 array, start with identity matrices
#     M = np.asarray([[[1,0],[0,1]] for i in range(lam.size)])
    
#     # gamma for incident material
#     gg0 = _gamma(n0,alpha0,pol)
    
#     # incidence angle for substrate material
#     alphas = np.arcsin(n0/ns*np.sin(alpha0))
    
#     # gamma for substrate material
#     ggs = _gamma(ns,alphas,pol)
    
#     # calculate the nx2x2 array for each layer, multiply the 2x2's, for each layer
#     for nd in ndpairs:
#         nj = nd[0]
#         dj = nd[1]
        
#         alphaj = np.arcsin(n0/nj*np.sin(alpha0))
#         mj = _mmj(lam,nj,dj,alphaj,pol)
        
#         M = np.matmul(mj,M)
        
#         mlist = list(M)

#     # calculate the (amplitude) reflection coeffecient given the total matrix
#     ref = [(g0*m[0,0]+gs*g0*m[0,1]-m[1,0]-gs*m[1,1])/(g0*m[0,0]+gs*g0*m[0,1]+m[1,0]+gs*m[1,1])\
#                for g0,gs,m in zip(gg0,ggs,mlist)]
    
#     # return the intensity reflectance
#     return ref

# display gamma correction
def gammacor(x):
    '''
    Gamma corrects the R, G, or B value

    Parameters
    ----------
    x : float
        R, G, or B value between 0 and 1.

    Returns
    -------
    float
        gamma-corrected R, G, or B value.

    '''
    return np.piecewise(x, [x <= 0.00313, x > 0.00313],
                  [lambda v: 12.92 * v,
                  lambda v: 1.055 * v ** (1 / 2.4) - 0.055])

#calculate the 
def refrgb(lam,incspec,n0,alpha0,pol,ns,*ndpairs):
    '''

    Parameters
    ----------
    
    lam : 1D numpy array of shape (nlam,) 
        wavelength range in nanometers.
    incspec : 1D numpy array of shape(nlam,)
        Incident power spectrum.
    n0 : 1D numpy array of shape (nlam,)
        refractive index of incident material at the values of lam
    alpha0 : int or float
        angle of incidence for incident layer, default normal incidence
    pol : string
        'TE', 'TM', or 'mixed' (50% each), otherwise normal incidence.
    ns : 1D numpy array of shape (nlam,)
        refractive index of substrate material at the values of lam
    *ndpairs : variable number of 2 element lists of [1d numpy array,int/float ]
       The RI n(λ) and thickness d (in nm) of each layer.

    Returns
    -------
    Numpy array of shape (3,)
        RGB values between 0 and 1.

    '''
           
    #calculate the intensity reflectance spectrum 
    spec = refspec(lam,n0,alpha0,pol,ns,*ndpairs)*incspec
    
    # 'integrate' / dot to get the CIE x, y, and z components
    x = sum(ciex*spec)
    y = sum(ciey*spec)
    z = sum(ciez*spec)
    
    # make into an array and divide (factor of 5 empirically agrees with paper)
    xyz = np.asarray([x,y,z])/5
    
    # convert xyz to rgb
    rgb = np.dot(xyz,xyz_to_rgb)
    
    # perform gamma correction
    rgb = gammacor(rgb)
    
    # clip values less than 0 or greater than 1
    rgb[rgb<0]=0
    rgb[rgb>1]=1

    return rgb    