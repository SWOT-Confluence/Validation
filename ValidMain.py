# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 16:14:14 2021

@author: coss.31
"""

from os import listdir
from os.path import isfile, join
import inputGncREADER as Rin
import validation as valid
import W2cdfvalid as ncO
import numpy as np

EMPTY=-9999

gagedir= "C:/Users/coss.31/OneDrive - The Ohio State University/Documents/PYfun/valid/gageQ"
SWOTqdir="C:/Users/coss.31/OneDrive - The Ohio State University/Documents/PYfun/valid/SWOTq/offline"

#Find validation data NC files--This will need updates when more agencies are added!!
onlyfiles = [f for f in listdir(gagedir) if isfile(join(gagedir, f))]
for file in onlyfiles:
    if 'USGS' in file:
        readnc=gagedir+'/'+file
        USGS_reach, USGSq,USGSt=Rin.read(readnc)#read gage file
    elif 'GRDC' in file:
            readnc=gagedir+'/'+file
            GRDC_reach, GRDCq,GRDCt=Rin.read(readnc)#read gage file       
# find input SWOT files taht require validation
onlyfiles = [f for f in listdir(SWOTqdir) if isfile(join(SWOTqdir, f))]

#loop over reach IDs in SWOT input cheching validation files for matches
validout=()            
for file in onlyfiles:
    #check input reach name against reaches in validation datasets
    IDstr=file[0:11]
         
    L1=IDstr in USGS_reach    
    L2=IDstr in GRDC_reach   
    if L1 or L2:
        SWOT_file=SWOTqdir+'/'+file
        SWt,SWq=Rin.SWOTread(SWOT_file)
        if L1:
            Target = np.where(USGS_reach==IDstr)  
            Vt=USGSt[Target,:]
            Vq=USGSq[Target,:]          
        elif L2:
            Target =  np.where(GRDC_reach==IDstr)  
            Vt=GRDCt[:,Target]
            Vt = Vt.astype(int)+366 # Vt in matlab tatetime +366 converts to Python ordinal
            Vq=GRDCq[:,Target]   
             
        Vout= valid.stats(SWt,SWq,Vt,Vq,IDstr)
        ncO.write(Vout,IDstr,SWt)
    else:
        ncO.write(EMPTY,IDstr,EMPTY)
       
    