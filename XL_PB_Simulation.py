#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 16:26:34 2018

@author: bryanquigley

X-ray Luminescence Simulation

"""

from __future__ import division

import sys
sys.path.insert(0,"/home/bquigley/Research/XFXL Simulation 1.0.0")

import numpy as np
import matplotlib.image
import Surface_Radiance
import os


if not os.path.exists('XL_physics_images'):
    os.makedirs('XL_physics_images')
if not os.path.exists('XL_Simulated_images'):
    os.makedirs('XL_Simulated_images')

#create capillary tube mask
tubemask = np.zeros((5000,512), dtype=np.float64)
voxelsize =  0.0100 # 0.01 mm
voxelsizex = 0.21586 #mm
tubediam = 0.59436 # 0.6 mm
spacing = 2.262 # 2.262 mm (2.862 mm if you measure from center-to-center of tubes)
tube1center = 255 - (spacing/2)/voxelsizex - (tubediam/2)/voxelsizex
tube2center = 255 + (spacing/2)/voxelsizex + (tubediam/2)/voxelsizex
Rsquared = (tubediam/2)**2

for i in xrange(5000):
    for j in xrange(512):
        if j < 255:

            if ((i-2500)*voxelsize) ** 2 + ((j-tube1center)*voxelsizex) ** 2 <= Rsquared:
                tubemask[i,j] = 1

        else:
            if ((i-2500)*voxelsize) ** 2 + ((j-tube2center)*voxelsizex) ** 2 <= Rsquared:
                tubemask[i,j] = 1

matplotlib.image.imsave('XL_physics_images/tubemask.jpg', tubemask)

#create gel phantom around tubes
phantommask = np.zeros((5000,512), dtype=np.float64)
for j in xrange(5000):
    if j>=1700:
        phantommask[j,:]=1

for i in xrange(5000):
    for j in xrange(512):
        if tubemask[i,j]==1:
            phantommask[i,j]=2

matplotlib.image.imsave('XL_physics_images/phantommask.jpg', phantommask)

xdim = 0.21586 # pixel size in object space 0.02297255 mm
ydim = 0.200 #pencilbeam step size in y direction mm
PBwidth = 0.156 # pencil beam width in mm. assume square shape
     
#attenuation map
AttMap = np.zeros((250,512), dtype=np.float64)
PBattenuation_water = np.zeros((512), dtype=np.float64)
PBattenuation_air = np.zeros((512), dtype=np.float64)
PBenergy = 17.4 # 17.4 keV
mac_water = 1.258568 # cm^2/g for 17.4 keV x-ray photons
rho_water = 1.000E+00 #g/cm^3 density of water(phantom)
mac_air = 1.212672 # cm^2/g for 17.4 keV x-ray photons
rho_air = 1.205E-03 #g/cm^3 density of air
mac_NPs = 77.66 # cm^2/g for Y2O3:Eu3+
rho_NPs = 1.00/(4*1.415*1.113366) # g/cm^3

meac_water = 0.978624 #cm^/g for 17.4 keV x-ray photons
meac_air = 0.952352 #cm^/g for 17.4 keV x-ray photons
meac_NP = 36.91 #cm^2/g for 17.4 keV x-rayt photons

density_y=np.zeros((250,512), dtype=np.float64)

for i in xrange(250):
    for j in xrange(512):
        density_y[i,j] = sum(tubemask[0+i*20:i*20+20, j])/20


muMap = np.zeros((250,512), dtype=np.float64)
meacMap = np.zeros((250,512), dtype=np.float64)
    
for y in xrange(250):
    for x in xrange(512):
        if y< 85:
            muMap[y,x] = mac_air
            meacMap[y,x] = meac_air
        else:
            muMap[y,x] = mac_NPs*density_y[y,x]*rho_NPs+mac_water*(1-density_y[y,x])*rho_water
            meacMap[y,x] = meac_NP*density_y[y,x]+meac_water*(1-density_y[y,x])
                
AttMap[:,0]=np.exp(-1*muMap[:,0]*(xdim*0.1-3.00443))    

for y in xrange(250):
    for x in xrange(511):
        AttMap[y,x+1]=AttMap[y,x]*np.exp(-1*muMap[y,x+1]*0.1*xdim)

I0 = 1.02728E+11 #Incident photon flux density Numbers/s/cm^2
t = 4 #dwell time 600 s
E = 17.4 # keV  

XrayDose = np.zeros((250,512), dtype=np.float64)

for i in xrange(250):
    for j in xrange(512):
        XrayDose[i,j] = I0 * E * t * meacMap[i,j] * AttMap[i,j]


matplotlib.image.imsave('XL_physics_images/XrayDose.jpg', XrayDose)

#X-ray Luminescence: Calculate [D_voxel * rho_phan * C_np * alpha_np]

rho_phan = 1 #g/cm^3
C_NP = 1/(4*1.415*1.113366) # g/cm^3 (estimate)
alpha_NP = 1.06 #photons/keV/g/cm^3

xrayluminescence = np.zeros((250,512), dtype=np.float64)
   
for i in xrange(250):
    for j in xrange(512):
        xrayluminescence[i,j] = XrayDose[i,j] * rho_phan * C_NP * alpha_NP * density_y[i,j] * xdim*0.1 * (PBwidth*.1)**2
       
matplotlib.image.imsave('XL_physics_images/xrayluminescence.jpg', xrayluminescence)

surfaceluminescence = np.zeros((512,512,250), dtype=np.float64)

for z in xrange(250):
    for i in xrange(512):
        if xrayluminescence[z,i] > 0:
            for j in xrange(512):
                for k in xrange(512):
                    r=xdim*0.1*np.sqrt((k-i)**2 + (j-255)**2)
                    surfaceluminescence[j,k,z]=surfaceluminescence[j,k,z]+Surface_Radiance.Nsurfrad(r, (z-85)*ydim*0.1, 0.25, 5)*xrayluminescence[z,i]
    print('Pencil beam position ' + str(z) + ' finished')                
            
surflum124=surfaceluminescence[:,:,124] 
matplotlib.image.imsave('XL_physics_images/SurfaceLuminescence124.jpg', surflum124)
 
#detected surface radiance:
mag = 16/215.86 #magnification
N = 2.8 # f/number
ps = 0.021586 # cm pixel dimension in image space
epsilon = 0.9

DetectedPhotons = np.zeros((512,512,250), dtype=np.float64)

for z in xrange(250):
    for i in xrange(512):
        for j in xrange(512):
            DetectedPhotons[i,j,z] = np.pi / 4 / ((1+mag)**2) / (N**2) * ps**2 * epsilon * surfaceluminescence[i,j,z]

RandomDetectorSignal=np.zeros((512,512,250), dtype=np.float64)

dc = 0.012 # e-/pixel/s dark current
M = 1000 # Electron Multipling gain
BE = 0.005 # e-/pixel EMCCD-amplified background events
readnoise = 0.1 # -e/pixel readnoise

#Random number distribution
for z in xrange(113,136):
    for i in xrange(512):
        for j in xrange(512):
            RandomDetectorSignal[i,j,z]=np.random.poisson(DetectedPhotons[i,j,z])*M + np.random.normal(0, BE)*M + np.random.normal(0, dc*t)*M + np.random.normal(0, readnoise)
    print('finished ' + str(z) + 'th itteration of randomized detector noise')
    np.save('XL_Simulated_images/Sim_XL_position_'+str(z-112),RandomDetectorSignal[:,:,z])
    matplotlib.image.imsave('XL_Simulated_images/Sim_XL_position_'+str(z-112)+'.tif', np.array(RandomDetectorSignal[:,:,z]), cmap='gray')
    
