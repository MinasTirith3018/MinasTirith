#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 18:12:36 2017

@author: ao2017
"""
cd G:\IAO/R
cd 2
import os
import time
import subprocess


fol=['I','V','R']
date=[20071025,20071027,20071028]
num=range(1,7)

j=0
k=20071025
i=2
for j in range(3):
for k in date:
    for i in num:
        ct=0
        count=0
        while (ct<3):
            os.chdir('/home/ao2017/toshiba/IAO/'+str(fol[j])+'/'+ (str(k))+'/'+(str(i)))
            ct+=1
        time.sleep(1)
    
        while (count<5):
                dirlist = os.listdir('.')
                dirlist = list(filter(lambda x: x[0:3] == '17P', dirlist))
                output = '\n'.join(dirlist)
                f = open('light.txt', 'w')
                f.write(output)
                count+=1

        time.sleep(1)                    
    
        subprocess.call(["python", "filter_2.py", "--list", "light.txt", "--output", "2_light_{:}.txt".format(i)])
        subprocess.call(["python", "cometregister.py", "--list", "light_{:}.txt".format(i), "--output", "offset_{:}.txt".format(i)])
        subprocess.call(["python", "preprocessor.py", "--masterdark", "masterdark_300.0_253.13_0.fits", "--masterflat", "masterflat_R_253.13_0.fits", "--light", "light_{:}.txt".format(i), "--fits_header_ccdtemp", "SET-TMP", "--offset", "offset_{:}.txt".format(i)])

i=1

#%%
fol=['I','V','R']
dat=[1025,1027,1028]


    for k in dat:
        ct=0
        count=0
        while (ct<3):
                os.chdir('/home/ao2017/toshiba/IAO/R/')
                ct+=1
        time.sleep(1)
        
        while (count<5):
                    dirlist = os.listdir('.')
                    dirlist = list(filter(lambda x: x[0:3] == str(k), dirlist))
                    output = '\n'.join(dirlist)
                    f = open('light_{:}.txt'.format(k), 'w')
                    f.write(output)
                    count+=1
                    
                    
        ls 1025*.fit > light_1025.txt            
        k=1025
        time.sleep(1)                    
        subprocess.call(["python", "cometregister.py", "--list", "light_{:}.txt".format(k), "--output", "offset_{:}.txt".format(k)])
        subprocess.call(["python", "preprocessor.py", "--light", "light_{:}.txt".format(k), "--fits_header_ccdtemp", "SET-TMP", "--offset", "offset_{:}.txt".format(k)])

#%% for measure the expose time
fol=['I','V','R']
date=[20071025,20071027,20071028]
num=range(1,7)


for j in range(3):
    for k in date:
        for i in num:
            ct=0
            count=0
            while (ct<3):
                os.chdir('/home/ao2017/toshiba/IAO/'+str(fol[j])+'/'+ (str(k))+'/'+(str(i)))
                ct+=1
            time.sleep(1)
        
            while (count<5):
                    dirlist = os.listdir('.')
                    dirlist = list(filter(lambda x: x[0:3] == '17P', dirlist))
                    output = '\n'.join(dirlist)
                    f = open('expos.txt', 'w')
                    f.write(output)
                    count+=1
    
            time.sleep(1)                    
            subprocess.call(["python", "preprocessor.py", "--light", "expos.txt", "--fits_header_ccdtemp", "SET-TMP"])
    
#%% Mesure the expose time 2

text_file = open("light.txt",'r')
lines = text_file.readlines()
len(lines)

exptime=0.
for i in lines:
            i=str(i).replace("\n","")
            img = CCDData.read(str(i), hdu=0, unit='u.electron/u.s')
            exptime+=img.header['EXPTIME']

exptime

#%% 다크 노출시간 체크 및 txt파일로 기록

cd ../20071028

text_file = open("dark.txt",'r')
lines = text_file.readlines()
len(lines)

darklist=[]
f=open("dark1s.txt","w")
for i in lines:
            i=str(i).replace("\n","")
            img = CCDData.read(str(i), hdu=0, unit='u.electron/u.s')
            if (img.header['EXPTIME']==1.):
                bias=str(i)
                f.write(bias+'\n')
                darklist.append(str(i))
                
f.close()
darklist

#%% 폴더별 노출시간 기록

cd G:\IAO/R

fol=['I','V','R']
date=[20071025,20071027,20071028]
num=range(1,7)

for k in date:
     for i in num:
        ct=0
        while (ct<3):
                os.chdir('G:\IAO/R/' + (str(k))+'/'+(str(i)))
                ct+=1
        time.sleep(1)
        while (ct<8):
            dirlist = os.listdir('.')
            dirlist = list(filter(lambda x: x[0:3] == '17P', dirlist))
            output = '\n'.join(dirlist)
            f = open('maxim_light.txt', 'w')
            f.write(output)
            ct+=1
 f.close()

expt=np.zeros((14,3))
count=0
for k in date:
    for j in range(1,7):
        exptime=0.
        os.chdir('G:\IAO/R/' + (str(k))+'/'+(str(j)))        
        time.sleep(1)        
        text_file= open('maxim_light.txt', 'r')
        lines=text_file.readlines()        
        time.sleep(1)  
        for i in lines:
                i=str(i).replace('\n',"")
                img=CCDData.read(str(i), hdu=0, unit='u.electron/u.s')
                exptime+=img.header['EXPTIME']
        expt[count,0]=k
        expt[count,1]=j
        expt[count,2]=exptime
        count+=1
        
cd ../
np.savetxt('t_exp_I(1).csv',expt[:,0])
np.savetxt('t_exp_I(2).csv',expt[:,1])
np.savetxt('t_exp_I.csv',expt[:,2])
