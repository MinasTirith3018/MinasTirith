# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 01:01:23 2017

@author: Canopus
"""
cd
cd Desktop/Maxim-3

from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from ccdproc import CCDData
from astropy.convolution import Gaussian2DKernel
from astropy.wcs import WCS
from photutils import detect_sources, source_properties, properties_table
import time
import subprocess
import os

 day=[1025,1027,1028]

        #y_cent= np.argwhere(img.data == img.data.max())[0][0]#원 그리는 센터
        #x_cent= np.argwhere(img.data == img.data.max())[0][1]
        #plt.plot(x_cent, y_cent, marker='x', ms=10, color='r')
   
#%%
  


cd ../Maxim-V
    
jd_v=np.zeros(14)
max_v=np.zeros(14)
snap_v=np.zeros(14)
area_v=np.zeros(14)
area2_v=np.zeros(14)
TotalMag_v=np.zeros(14)
TotalMag2_v=np.zeros(14)
exp_v=np.zeros(14)
utc_v=[]
frame_v=[]


       
        
exp_v=np.array([281.,120.,150,50,100,94,870,940,600,400,200,255,920,360]) #real exp by maxim

i=1
j=1025
count=0
for j in day:
    for i in range(1,7):
        img = CCDData.read('V-{:}-{:}(2).fit'.format(j,i), hdu=0, unit='u.electron/u.s')
        exp_v[count]=img.header['EXPOSURE']
        pic=img.data/exp_v[count]
        time.sleep(1)
        jd_v[count]=img.header['JD']
        max_v[count]=pic.max()
        snap_v[count]=img.header['SNAPSHOT']
        utc_v.append(img.header['UT'])
        frame_v.append(img.header['FRAMEID'])
        lev=np.linspace(1+pic.min(),pic.max(),31)
        lev2=np.linspace(-2.5*log10(pic.max()),-2.5*log10(1.+pic.min()),31)
        midpoint1=(min(lev2)+2*max(lev2))/3
        #midpoint2=(min(lev2))/3
        
        comet=np.zeros(np.shape(pic))
        #comet2=np.zeros(np.shape(pic))
        area=0.
        #area2=0.
        for x in range(1024):
            for y in range(1074):
                if (pic[y,x] >0):
                    if (-2.5*log10(pic[y,x]) <= midpoint1) : 
                        comet[y,x]=pic[y,x]
                        area+=1.
                        TotalMag_v[count]+=comet[y,x]
        time.sleep(1)
        area_v[count]=area      
        #for x in range(1024):
        #    for y in range(1074):
        #        if (pic[y,x] >0):
        #            if (-2.5*log10(pic[y,x]) <= midpoint2) : 
        #                comet2[y,x]=pic[y,x]
        #                area2+=1.
        #                TotalMag2_v[count]+=comet2[y,x]
        #area2_v[count]=area2
        #plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(comet2),lev2),plt.colorbar(),plt.axis('equal'),plt.title('Coma : V-{:}-{:}'.format(j,i)),plt.savefig('V-{:}-{:} (2)coma.png'.format(j,i)) 
        plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(comet),lev2),plt.colorbar(),plt.axis('equal'),plt.title('Coma : V-{:}-{:}'.format(j,i)),plt.savefig('V-{:}-{:} (1)coma.png'.format(j,i))      
        plt.figure(figsize=(6,5)),plt.contourf(pic,lev),plt.colorbar(),plt.axis('equal'),plt.title('V-{:}-{:}'.format(j,i)),plt.savefig('V-{:}-{:} (1).png'.format(j,i))
        plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(pic),lev2),plt.colorbar(),plt.axis('equal'),plt.title('V-{:}-{:}'.format(j,i)),plt.savefig('V-{:}-{:} (2).png'.format(j,i))
        count+=1



       
mag_v=-2.5*log10(max_v)
plt.plot(jd_v,max_v)
plt.figure(figsize=(8,6)),plt.plot(jd_v,mag_v, 'b+-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log( intensity'),plt.title('V band-Mag'),plt.savefig('V-Mag.png')
plt.figure(figsize=(8,6)),plt.plot(jd_v,0.72*sqrt(area_v/pi),'o-'),plt.xlabel('Julian Day'),plt.ylabel('Radius/arcsec'),plt.title('V-band Coma radius'),plt.savefig('V-Radius.png')
plt.figure(figsize=(8,6)),plt.plot(jd_v,0.72*sqrt(area2_v/pi),'o-'),plt.xlabel('Julian Day'),plt.ylabel('Radius/arcsec'),plt.title('V-band Coma radius'),plt.savefig('V-Radius-2.png')

meanMag_v=-2.5*log10(TotalMag_v/(area_v))
TMag_v=-2.5*log10(TotalMag_v)
plt.figure(figsize=(8,6)),plt.plot(jd_v,meanMag_v,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Intensity/pix)'),plt.title('V band Mean-mag'),plt.savefig('V-MeanMag.png')
plt.figure(figsize=(8,6)),plt.plot(jd_v,TMag_v,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Total Intensity)'),plt.title('V band Total Mag'),plt.savefig('V-TotalMag.png')

meanMag2_v=-2.5*log10(TotalMag2_v/(area2_v))
TMag2_v=-2.5*log10(TotalMag2_v)
plt.figure(figsize=(8,6)),plt.plot(jd_v,meanMag2_v,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Intensity/pix)'),plt.title('V band Mean-mag'),plt.savefig('V-MeanMag-2.png')
plt.figure(figsize=(8,6)),plt.plot(jd_v,TMag2_v,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Total Intensity)'),plt.title('V band Total Mag'),plt.savefig('V-TotalMag-2.png')


np.savetxt('JulianDay.csv',jd_v)
np.savetxt('area.csv',area_v,fmt='%.2f')
#np.savetxt('area2.csv',area2_v,fmt='%.2f')
np.savetxt('radius-pix.csv',sqrt(area_v/pi))
#np.savetxt('radius-pix2.csv',sqrt(area2_v/pi))
np.savetxt('maxvalue.csv',max_v)
np.savetxt('maxvalue magnitude.csv',mag_v)
np.savetxt('Total intensity.csv',TotalMag_v)
np.savetxt('Mean Coma magnitude.csv',meanMag_v)
np.savetxt('Total Coma magnitude.csv', TMag_v)
#np.savetxt('Total intensity2.csv',TotalMag2_v,fmt='%.2f')
#np.savetxt('Mean Coma magnitude2.csv',meanMag2_v)
#np.savetxt('Total Coma magnitude2.csv', TMag2_v)
np.savetxt('snapnumber.csv',snap_v)
np.savetxt('t_exp_V.csv',exp_v)

output= open('UTC.csv','w')
for i in utc_v: output.write(i+'\n')
output.close()

output= open('Reference.csv','w')
for i in frame_v: output.write(i+'\n')
output.close()


#%%

cd ../../Maxim-3/Maxim-I  

    day=[1025,1027,1028]


    
jd_i=np.zeros(14)
max_i=np.zeros(14)
snap_i=np.zeros(14)
area_i=np.zeros(14)
area2_i=np.zeros(14)
TotalMag_i=np.zeros(14)
TotalMag2_i=np.zeros(14)
exp_i=np.zeros(14)
utc_i=[]
frame_i=[]


       
        
exp_i=np.array([281.,120.,150,50,100,94,870,940,600,400,200,255,920,360]) #real exp by maxim
exp_i

count=0
for j in day:
    for i in range(1,7):
        img = CCDData.read('I-{:}-{:}(2).fit'.format(j,i), hdu=0, unit='u.electron/u.s')
        exp_i[count]=img.header['EXPOSURE']
        pic=img.data/exp_i[count]
        time.sleep(1)
        jd_i[count]=img.header['JD']
        max_i[count]=pic.max()
        snap_i[count]=img.header['SNAPSHOT']
        utc_i.append(img.header['UT'])
        frame_i.append(img.header['FRAMEID'])
        
        lev=np.linspace(1.+pic.min(),pic.max(),31)
        lev2=np.linspace(-2.5*log10(pic.max()),-2.5*log10(1.+pic.min()),31)
        midpoint1=(min(lev2)+2*max(lev2))/3
        midpoint2=(min(lev2))/3
        
        comet=np.zeros(np.shape(pic))
        comet2=np.zeros(np.shape(pic))
        area=0.
        area2=0.
        for x in range(1024):
            for y in range(1074):
                if (pic[y,x] >0):
                    if (-2.5*log10(pic[y,x]) <= midpoint1) : 
                        comet[y,x]=pic[y,x]
                        area+=1.
                        TotalMag_i[count]+=comet[y,x]
        time.sleep(1)
        area_i[count]=area
        #for x in range(1024):
        #    for y in range(1074):
        #        if (pic[y,x] >0):
        #            if (-2.5*log10(pic[y,x]) <= midpoint2) : 
        #                comet2[y,x]=pic[y,x]
        #                area2+=1.
        #                TotalMag2_i[count]+=comet2[y,x]
        #area2_i[count]=area2
        
        plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(comet),lev2),plt.colorbar(),plt.axis('equal'),plt.title('Coma : I-{:}-{:}'.format(j,i)),plt.savefig('I-{:}-{:} (1)coma.png'.format(j,i))      
        #plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(comet2),lev2),plt.colorbar(),plt.axis('equal'),plt.title('Coma : I-{:}-{:}'.format(j,i)),plt.savefig('I-{:}-{:} (2)coma.png'.format(j,i)) 
        plt.figure(figsize=(6,5)),plt.contourf(pic,lev),plt.colorbar(),plt.axis('equal'),plt.title('I-{:}-{:}'.format(j,i)),plt.savefig('I-{:}-{:} (1).png'.format(j,i))
        plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(pic),lev2),plt.colorbar(),plt.axis('equal'),plt.title('I-{:}-{:}'.format(j,i)),plt.savefig('I-{:}-{:} (2).png'.format(j,i))
        count+=1



mag_i=-2.5*log10(max_i)
plt.plot(jd_i,max_i)
plt.figure(figsize=(8,6)),plt.plot(jd_i,mag_i, 'b+-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log( intensity'),plt.title('I band-Mag'),plt.savefig('I-Mag.png')
plt.figure(figsize=(8,6)),plt.plot(jd_i,0.72*sqrt(area_i/pi),'o-'),plt.xlabel('Julian Day'),plt.ylabel('Radius/arcsec'),plt.title('I-band Coma radius'),plt.savefig('I-Radius.png')
#plt.figure(figsize=(8,6)),plt.plot(jd_i,0.72*sqrt(area2_i/pi),'o-'),plt.xlabel('Julian Day'),plt.ylabel('Radius/arcsec'),plt.title('I-band Coma radius'),plt.savefig('I-Radius-2.png')

meanMag_i=-2.5*log10(TotalMag_i/(area_i))
TMag_i=-2.5*log10(TotalMag_i)
plt.figure(figsize=(8,6)),plt.plot(jd_i,meanMag_i,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Intensity/pix)'),plt.title('I band Mean-mag'),plt.savefig('I-MeanMag.png')
plt.figure(figsize=(8,6)),plt.plot(jd_i,TMag_i,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Total Intensity)'),plt.title('I band Total Mag'),plt.savefig('I-TotalMag.png')

#meanMag2_i=-2.5*log10(TotalMag2_i/(area2_i))
#TMag2_i=-2.5*log10(TotalMag2_i)
#plt.figure(figsize=(8,6)),plt.plot(jd_i,meanMag2_i,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Intensity/pix)'),plt.title('I band Mean-mag'),plt.savefig('I-MeanMag-2.png')
#plt.figure(figsize=(8,6)),plt.plot(jd_i,TMag2_i,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Total Intensity)'),plt.title('I band Total Mag'),plt.savefig('I-TotalMag-2.png')


np.savetxt('JulianDay.csv',jd_i)
np.savetxt('area.csv',area_i,fmt='%.2f')
#np.savetxt('area2.csv',area2_i,fmt='%.2f')
np.savetxt('radius-pix.csv',sqrt(area_i/pi))
#np.savetxt('radius-pix2.csv',sqrt(area2_i/pi))
np.savetxt('maxvalue.csv',max_i)
np.savetxt('maxvalue magnitude.csv',mag_i)
np.savetxt('Total intensity.csv',TotalMag_i)
np.savetxt('Mean Coma magnitude.csv',meanMag_i)
np.savetxt('Total Coma magnitude.csv', TMag_i)
#np.savetxt('Total intensity2.csv',TotalMag2_i,fmt='%.2f')
#np.savetxt('Mean Coma magnitude2.csv',meanMag2_i)
#np.savetxt('Total Coma magnitude2.csv', TMag2_i)
np.savetxt('snapnumber.csv',snap_i)
np.savetxt('t_exp_I.csv',exp_i)

output= open('UTC.csv','w')
for i in utc_i: output.write(i+'\n')
output.close()

output= open('Reference.csv','w')
for i in frame_i: output.write(i+'\n')
output.close()







#%%

cd ../../Maxim-3/Maxim-R  

    day=[1025,1027,1028]


    
jd_r=np.zeros(14)
max_r=np.zeros(14)
snap_r=np.zeros(14)
area_r=np.zeros(14)
area2_r=np.zeros(14)
TotalMag_r=np.zeros(14)
TotalMag2_r=np.zeros(14)
exp_r=np.zeros(14)
utc_r=[]
frame_r=[]


       
        
exp_r=np.array([281.,120.,150,50,100,94,870,940,600,400,200,255,920,360]) #real exp by maxim
exp_r 

count=0
for j in day:
    for i in range(1,7):
        img = CCDData.read('R-{:}-{:}(2).fit'.format(j,i), hdu=0, unit='u.electron/u.s')
        exp_r[count]=img.header['EXPOSURE']
        pic=img.data/exp_r[count]
        time.sleep(1)
        jd_r[count]=img.header['JD']
        max_r[count]=pic.max()
        snap_r[count]=img.header['SNAPSHOT']
        utc_r.append(img.header['UT'])
        frame_r.append(img.header['FRAMEID'])

        lev=np.linspace(1+pic.min(),pic.max(),31)
        lev2=np.linspace(-2.5*log10(pic.max()),-2.5*log10(1.+pic.min()),31)
        midpoint1=(min(lev2)+2*max(lev2))/3
        midpoint2=(min(lev2))/3
        
        comet=np.zeros(np.shape(pic))
        comet2=np.zeros(np.shape(pic))
        area=0.
        area2=0.
        for x in range(1024):
            for y in range(1074):
                if (pic[y,x] >0):
                    if (-2.5*log10(pic[y,x]) <= midpoint1) : 
                        comet[y,x]=pic[y,x]
                        area+=1.
                        TotalMag_r[count]+=comet[y,x]
        time.sleep(1)
        area_r[count]=area
       # for x in range(1024):
       #     for y in range(1074):
       #         if (pic[y,x] >0):
       #             if (-2.5*log10(pic[y,x]) <= midpoint2) : 
       #                 comet2[y,x]=pic[y,x]
       #                 area2+=1.
       #                 TotalMag2_r[count]+=comet2[y,x]
       # area2_r[count]=area2
        
        plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(comet),lev2),plt.colorbar(),plt.axis('equal'),plt.title('Coma : R-{:}-{:}'.format(j,i)),plt.savefig('R-{:}-{:} (1)coma.png'.format(j,i))      
       # plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(comet2),lev2),plt.colorbar(),plt.axis('equal'),plt.title('Coma : R-{:}-{:}'.format(j,i)),plt.savefig('R-{:}-{:} (2)coma.png'.format(j,i)) 
        plt.figure(figsize=(6,5)),plt.contourf(pic,lev),plt.colorbar(),plt.axis('equal'),plt.title('R-{:}-{:}'.format(j,i)),plt.savefig('R-{:}-{:} (1).png'.format(j,i))
        plt.figure(figsize=(6,5)),plt.contourf(-2.5*log10(pic),lev2),plt.colorbar(),plt.axis('equal'),plt.title('R-{:}-{:}'.format(j,i)),plt.savefig('R-{:}-{:} (2).png'.format(j,i))
        count+=1



mag_r=-2.5*log10(max_r)
plt.plot(jd_r,max_r)
plt.figure(figsize=(8,6)),plt.plot(jd_r,mag_r, 'b+-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log( intensity'),plt.title('R band-Mag'),plt.savefig('R-Mag.png')
plt.figure(figsize=(8,6)),plt.plot(jd_r,0.72*sqrt(area_r/pi),'o-'),plt.xlabel('Julian Day'),plt.ylabel('Radius/arcsec'),plt.title('R-band Coma radius'),plt.savefig('R-Radius.png')
#plt.figure(figsize=(8,6)),plt.plot(jd_r,0.72*sqrt(area2_r/pi),'o-'),plt.xlabel('Julian Day'),plt.ylabel('Radius/arcsec'),plt.title('R-band Coma radius'),plt.savefig('R-Radius-2.png')

meanMag_r=-2.5*log10(TotalMag_i/(area_r))
TMag_r=-2.5*log10(TotalMag_r)
plt.figure(figsize=(8,6)),plt.plot(jd_r,meanMag_r,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Intensity/pix)'),plt.title('R band Mean-mag'),plt.savefig('R-MeanMag.png')
plt.figure(figsize=(8,6)),plt.plot(jd_r,TMag_r,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Total Intensity)'),plt.title('R band Total Mag'),plt.savefig('R-TotalMag.png')

#meanMag2_r=-2.5*log10(TotalMag2_r/(area2_r))
#TMag2_r=-2.5*log10(TotalMag2_r)
#plt.figure(figsize=(8,6)),plt.plot(jd_r,meanMag2_r,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Intensity/pix)'),plt.title('R band Mean-mag'),plt.savefig('R-MeanMag-2.png')
#plt.figure(figsize=(8,6)),plt.plot(jd_r,TMag2_r,'o-'),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(Total Intensity)'),plt.title('R band Total Mag'),plt.savefig('R-TotalMag-2.png')


np.savetxt('JulianDay.csv',jd_r)
np.savetxt('area.csv',area_r,fmt='%.2f')
#np.savetxt('area2.csv',area2_r,fmt='%.2f')
np.savetxt('radius-pix.csv',sqrt(area_r/pi))
#np.savetxt('radius-pix2.csv',sqrt(area2_r/pi))
np.savetxt('maxvalue.csv',max_r,)
np.savetxt('maxvalue magnitude.csv',mag_r)
np.savetxt('Total intensity.csv',TotalMag_r)
np.savetxt('Mean Coma magnitude.csv',meanMag_r)
np.savetxt('Total Coma magnitude.csv', TMag_r)
#np.savetxt('Total intensity2.csv',TotalMag2_r,fmt)
#np.savetxt('Mean Coma magnitude2.csv',meanMag2_r)
#np.savetxt('Total Coma magnitude2.csv', TMag2_r)
np.savetxt('snapnumber.csv',snap_r)
np.savetxt('t_exp_R.csv',exp_r)

output= open('UTC.csv','w')
for i in utc_r: output.write(i+'\n')
output.close()

output= open('Reference.csv','w')
for i in frame_r: output.write(i+'\n')
output.close()


#%% linear fitting 및 상관계수 그리는 영역(최종 결과물에는 안 씀)


cd ../
plt.figure(figsize=(8,6)),plt.plot(jd_v,mag_v, 'y+-',label='V band'),plt.plot(jd_r,mag_r, 'r+-',label='R band'),plt.plot(jd_i,mag_i, 'k+-',label='I band'),plt.legend(),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(intensity)'),plt.title('Apprent Mag'),plt.savefig('Mag2.png')
plt.figure(figsize=(8,6)),plt.plot(jd_v,meanMag_v, 'y+-',label='V band'),plt.plot(jd_r,meanMag_r, 'r+-',label='R band'),plt.plot(jd_i,meanMag_i, 'k+-',label='I band'),plt.legend(),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(intensity/pixel)'),plt.title('Mean Mag'),plt.savefig('Mag2_mean.png')
plt.figure(figsize=(8,6)),plt.plot(jd_v,TMag_v, 'y+-',label='V band'),plt.plot(jd_r,TMag_r, 'r+-',label='R band'),plt.plot(jd_i,TMag_i, 'k+-',label='I band'),plt.legend(),plt.xlabel('Julian Day'),plt.ylabel('-2.5log(sum intensity)'),plt.title('Sum Mag'),plt.savefig('Mag2_sum.png')



coeff_r, res_r, _, _, _ = np.polyfit(jd_r,0.72*sqrt(area_r/pi),1,full=True)
coeff_v, res_v, _, _, _ = np.polyfit(jd_v,0.72*sqrt(area_v/pi),1,full=True)
coeff_i, res_i, _, _, _ = np.polyfit(jd_i,0.72*sqrt(area_i/pi),1,full=True)



plt.figure(figsize=(8,6)),plt.plot(jd_r,0.72*sqrt(area_r/pi),'ro',label='R band'),plt.plot(jd_v,0.72*sqrt(area_v/pi),'yo',label='V band'),plt.plot(jd_i,0.72*sqrt(area_i/pi),'ko',label='I band'),\
          plt.plot(jd_r,coeff_r[0]*jd_r+coeff_r[1],'r-'),plt.plot(jd_i,coeff_i[0]*jd_i+coeff_i[1],'k-'),plt.plot(jd_v,coeff_v[0]*jd_v+coeff_v[1],'y-'),\
          plt.legend(),plt.xlabel('Julian Day'),plt.ylabel('Radius/arcsec'),plt.title('Coma radius'),plt.savefig('Radius2.png')

np.savetxt('V_radius_coeff.txt',coeff_v)
np.savetxt('V_radius_res.txt',res_v)
np.savetxt('R_radius_coeff.txt',coeff_r)
np.savetxt('R_radius_res.txt',res_r)
np.savetxt('I_radius_coeff.txt',coeff_i)
np.savetxt('I_radius_res.txt',res_i)




coeff1_r, res1_r,_,_,_ = np.polyfit(jd_r,TMag_r,1,full=True)
coeff1_v, res1_v,_,_,_ = np.polyfit(jd_v,TMag_v,1,full=True)
coeff1_i, res1_i,_,_,_ = np.polyfit(jd_i,TMag_i,1,full=True)
plt.figure(figsize=(8,6)),plt.plot(jd_r,TMag_r,'ro',label='R band'),plt.plot(jd_v,TMag_v,'yo',label='V band'),plt.plot(jd_i,TMag_i,'ko',label='I band'),\
          plt.plot(jd_r,coeff1_r[0]*jd_r+coeff1_r[1],'r-'),plt.plot(jd_i,coeff1_i[0]*jd_i+coeff1_i[1],'k-'),plt.plot(jd_v,coeff1_v[0]*jd_v+coeff1_v[1],'y-'),\
          plt.legend(),plt.xlabel('Julian Day'),plt.ylabel('Magnitude'),plt.title('Total Magnitude'),plt.savefig('Total Mag.png')

np.savetxt('V_radius_coeff1.txt',coeff1_v)
np.savetxt('V_radius_res1.txt',res1_v)
np.savetxt('R_radius_coeff1.txt',coeff1_r)
np.savetxt('R_radius_res1.txt',res1_r)
np.savetxt('I_radius_coeff1.txt',coeff1_i)
np.savetxt('I_radius_res1.txt',res1_i)




coeff2_r, res2_r,_,_,_ = np.polyfit(jd_r,meanMag_r,1,full=True)
coeff2_v, res2_v,_,_,_ = np.polyfit(jd_v,meanMag_v,1,full=True)
coeff2_i, res2_i,_,_,_ = np.polyfit(jd_i,meanMag_i,1,full=True)
plt.figure(figsize=(8,6)),plt.plot(jd_r,meanMag_r,'ro',label='R band'),plt.plot(jd_v,meanMag_v,'yo',label='V band'),plt.plot(jd_i,meanMag_i,'ko',label='I band'),\
          plt.plot(jd_r,coeff2_r[0]*jd_r+coeff2_r[1],'r-'),plt.plot(jd_i,coeff2_i[0]*jd_i+coeff2_i[1],'k-'),plt.plot(jd_v,coeff2_v[0]*jd_v+coeff2_v[1],'y-'),\
          plt.legend(),plt.xlabel('Julian Day'),plt.ylabel('Magnitude'),plt.title('Mean pixel Magnitude'),plt.savefig('Mean Mag.png')

np.savetxt('V_radius_coeff2.txt',coeff2_v)
np.savetxt('V_radius_res2.txt',res2_v)
np.savetxt('R_radius_coeff2.txt',coeff2_r)
np.savetxt('R_radius_res2.txt',res2_r)
np.savetxt('I_radius_coeff2.txt',coeff2_i)
np.savetxt('I_radius_res2.txt',res2_i)

        