#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:49:52 2017

@author: ao2017
"""

from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from ccdproc import CCDData
from astropy.convolution import Gaussian2DKernel
from astropy.wcs import WCS
from photutils import detect_sources, source_properties, properties_table


from astropy.coordinates import SkyCoord
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from photutils import EllipticalAnnulus as EllAn
from photutils import CircularAnnulus as Circul
from photutils import aperture_photometry as APPHOT


img = CCDData.read('V-1027-1.fit', hdu=0, unit='u.electron/u.s')

np.argwhere(img == img.data.max())
np.shape(img.data)
lev=np.linspace(1,img.data.max(),50)
lev2=log10(lev)
lev2[0]=0.
plt.contourf(log10(img.data),lev2),plt.colorbar(),plt.axis('equal')


#x_sdss, y_sdss = radec_sdss.to_pixel(wcs)#적경과 적위를 x,y픽셀좌표로 바꿔줌.
x_cent= np.argwhere(img.data == img.data.max())[0][0]#원 그리는 센터
y_cent= np.argwhere(img.data == img.data.max())[0][1]

#Plot
fig, ax = plt.subplots(1,1, figsize=(12,12))
ax.imshow(img.data, vmin=0, vmax=img.data.max(), origin='lower')
ax.plot(y_cent, x_cent, marker='x', ms=20, color='r')


day=[1025,1027,1028]
num=range(1,7)
#%%
cd Maxim-V
exp_t=np.loadtxt('t_exp_V.csv')#노출시간을 불러냄



ct=0


for j in day:
    for i in num:
        image = CCDData.read('V-{:}-{:}(2).fit'.format(j,i), hdu=0, unit='u.electron/u.s')
        rad= 450
        img=image.data/exp_t[ct]
        x_cent= np.argwhere(img == img.max())[0][0]#원 그리는 센터
        y_cent= np.argwhere(img == img.max())[0][1]

        annul = [Circul(positions=(y_cent,x_cent), r_in=b0, r_out=b0+1) for b0 in range(1,rad+1)]
        phot = APPHOT(img, annul)
        time.sleep(1)            
        
        #Plot
    
        semimaj = np.arange(1, rad+1)
        counts = np.zeros(rad)
        dcounts = np.zeros(rad)
        for k in range(0, rad):
            count = phot[phot.colnames[k+3]]
            counts[k] = count
            # phot.colnames = column names = "aperture_sum_X"         
            dcount = count /(annul[k].area())# count per pixel   
            if(dcount >0):
                dcounts[k]=dcount
            else :
                dcounts[k]=1.e-1                  
        N_surf=dcounts/dcounts[0]
        plt.figure(figsize=(6,4)),plt.plot(semimaj, N_surf, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Normalized Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.xlim(1,1e+3),plt.ylim(1.e-4,1),plt.savefig('V-{:}-{:}-Normalize curve(2).png'.format(j,i))
        plt.figure(figsize=(6,6)),plt.plot(semimaj, N_surf, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Normalized Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.xlim(1,1e+3),plt.ylim(1.e-4,1),plt.savefig('V-{:}-{:}-Normalize equal curve(2).png'.format(j,i))
        plt.figure(figsize=(6,4)),plt.plot(semimaj, dcounts, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.savefig('V-{:}-{:}-Surface curve(2).png'.format(j,i))
        plt.figure(figsize=(6,6)),plt.plot(semimaj, dcounts, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.axis('equal'),plt.savefig('V-{:}-{:}-Surface equal curve(2).png'.format(j,i))
        time.sleep(2)
        np.savetxt('V-{:}-{:}-radial pixel.csv'.format(j,i),semimaj)
        np.savetxt('V-{:}-{:}-Surface Brightness.csv'.format(j,i),dcounts)
        np.savetxt('V-{:}-{:}-Normalize Brightness.csv'.format(j,i),N_surf)
        ct+=1


        fig, ax = plt.subplots(1,1, figsize=(8,8))
        ax.imshow(img, vmin=0, vmax=img.max(), origin='lower')
        [annul[k].plot(ax=ax, color='red') for k in arange(1,rad+1,20)],plt.savefig('V-{:}-{:}-Annul circle.png'.format(j,i))

        np.savetxt('V-{:}-{:}-counts.csv'.format(j,i-1),counts)
        plt.plot(semimaj,counts)
#%%

cd Maxim-I
exp_t=np.loadtxt('t_exp_I.csv')

exp_t

ct=0

for j in day:
    for i in num:
        image = CCDData.read('I-{:}-{:}(2).fit'.format(j,i), hdu=0, unit='u.electron/u.s')
        rad= 450
        img=image.data/exp_t[ct]
        x_cent= np.argwhere(img == img.max())[0][0]#원 그리는 센터
        y_cent= np.argwhere(img == img.max())[0][1]

        annul = [Circul(positions=(y_cent,x_cent), r_in=b0, r_out=b0+1) for b0 in range(1,rad+1)]
        phot = APPHOT(img, annul)
        time.sleep(1)            
        
        #Plot
    
        semimaj = np.arange(1, rad+1)
        counts = np.zeros(rad)
        dcounts = np.zeros(rad)
        for k in range(0, rad):
            count = phot[phot.colnames[k+3]]
            counts[k] = count
            # phot.colnames = column names = "aperture_sum_X"         
            dcount = count /(annul[k].area())# count per pixel             
            dcounts[k]=dcount
                   
        N_surf=dcounts/dcounts[0]
        plt.figure(figsize=(6,4)),plt.plot(semimaj, N_surf, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Normalized Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.xlim(1,1e+3),plt.ylim(1.e-4,1),plt.savefig('I-{:}-{:}-Normalize curve.png'.format(j,i))
        plt.figure(figsize=(6,6)),plt.plot(semimaj, N_surf, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Normalized Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.xlim(1,1e+3),plt.ylim(1.e-4,1),plt.savefig('I-{:}-{:}-Normalize equal curve.png'.format(j,i))
        plt.figure(figsize=(6,4)),plt.plot(semimaj, dcounts, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.savefig('I-{:}-{:}-Surface curve.png'.format(j,i))
        plt.figure(figsize=(6,6)),plt.plot(semimaj, dcounts, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.axis('equal'),plt.savefig('I-{:}-{:}-Surface equal curve.png'.format(j,i))
        time.sleep(2)
        np.savetxt('I-{:}-{:}-radial pixel.csv'.format(j,i),semimaj)
        np.savetxt('I-{:}-{:}-Surface Brightness.csv'.format(j,i),dcounts)
        np.savetxt('I-{:}-{:}-Normalize Brightness.csv'.format(j,i),N_surf)
        ct+=1


        fig, ax = plt.subplots(1,1, figsize=(8,8))
        ax.imshow(img, vmin=0, vmax=img.max(), origin='lower')
        [annul[k].plot(ax=ax, color='red') for k in arange(1,rad+1,20)],plt.savefig('I-{:}-{:}-Annul circle.png'.format(j,i))

        np.savetxt('I-{:}-{:}-counts.csv'.format(j,i-1),counts)
        plt.plot(semimaj,counts)


#%%
cd ../Maxim-R

exp_t=np.loadtxt('t_exp_R.csv')

exp_t

ct=0

for j in day:
    for i in num:
        image = CCDData.read('R-{:}-{:}(2).fit'.format(j,i), hdu=0, unit='u.electron/u.s')
        rad= 450
        img=image.data/exp_t[ct]
        x_cent= np.argwhere(img == img.max())[0][0]#원 그리는 센터
        y_cent= np.argwhere(img == img.max())[0][1]

        annul = [Circul(positions=(y_cent,x_cent), r_in=b0, r_out=b0+1) for b0 in range(1,rad+1)]
        phot = APPHOT(img, annul)
        time.sleep(1)            
        
        #Plot
    
        semimaj = np.arange(1, rad+1)
        counts = np.zeros(rad)
        dcounts = np.zeros(rad)
        for k in range(0, rad):
            count = phot[phot.colnames[k+3]]
            counts[k] = count
            # phot.colnames = column names = "aperture_sum_X"         
            dcount = count /(annul[k].area())# count per pixel             
            if(dcount >0):
                dcounts[k]=dcount
            else :
                dcounts[k]=1.e-2
                   
        N_surf=dcounts/dcounts[0]

        plt.figure(figsize=(6,4)),plt.plot(semimaj, N_surf, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Normalized Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.ylim(1.e-4,2),plt.savefig('R-{:}-{:}-Normalize curve.png'.format(j,i))
        plt.figure(figsize=(6,6)),plt.plot(semimaj, N_surf, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Normalized Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.xlim(1,1e+3),plt.ylim(1.e-4,2),plt.savefig('R-{:}-{:}-Normalize equal curve.png'.format(j,i))
        plt.figure(figsize=(6,4)),plt.plot(semimaj, dcounts, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.savefig('R-{:}-{:}-Surface curve.png'.format(j,i))
        plt.figure(figsize=(6,6)),plt.plot(semimaj, dcounts, ls=':', marker='x'),plt.xscale('log'),plt.yscale('log'),plt.ylabel('Surface Bright(ADU/pixel)'),plt.xlabel('log Radial axis (pixel)'),plt.grid(ls=":"),plt.axis('equal'),plt.savefig('R-{:}-{:}-Surface equal curve.png'.format(j,i))
        time.sleep(2)
       
        np.savetxt('R-{:}-{:}-radial pixel.csv'.format(j,i),semimaj)
        np.savetxt('R-{:}-{:}-Surface Brightness.csv'.format(j,i),dcounts)
        np.savetxt('R-{:}-{:}-Normalize Brightness.csv'.format(j,i),N_surf)
        ct+=1


        fig, ax = plt.subplots(1,1, figsize=(8,8))
        ax.imshow(img.data, vmin=0, vmax=img.data.max(), origin='lower')
        [annul[k].plot(ax=ax, color='red') for k in arange(1,rad+1,20)],plt.savefig('R-{:}-{:}-Annul circle.png'.format(j,i))

        np.savetxt('R-{:}-{:}-counts.csv'.format(j,i-1),counts)
        plt.plot(semimaj,counts)


