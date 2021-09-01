# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 09:30:08 2021

@author: coss.31
"""
import pandas as pd
import numpy as np
def read(IDXnc):
    from netCDF4 import Dataset
   
    ncf = Dataset(IDXnc)
    if 'USGS'in IDXnc:
        # need to pull Q, T, and reach identifiers
        #reach
        USGS_reach=np.ma.getdata(ncf.variables['reach_ID'][:])
        USGS_reach=["%.f" % number for number in USGS_reach]
        #Q
        USGSq = np.ma.getdata(ncf.variables['USGS_Q'][:])    
        #Qt
        USGSt = np.ma.getdata(ncf.variables['USGS_Qt'][:]) 
        return  USGS_reach, USGSq,  USGSt
    elif 'GRDC' in IDXnc:
         # need to pull Q, T, and reach identifiers
        #reach
        GRDC_reach=np.ma.getdata(ncf.variables['Reach_ID'][:])
        GRDC_reach=["%.f" % number for number in GRDC_reach]
        #Q
        GRDCq = np.ma.getdata(ncf.variables['GRDC_Q'][:])    
        #Qt
        GRDCt = np.ma.getdata(ncf.variables['GRDC_Qt'][:]) 
        
            
        return GRDC_reach, GRDCq,  GRDCt
def SWOTread(IDXnc):
    from os import listdir
    from netCDF4 import Dataset
    from os.path import isfile, join
    from datetime import datetime
     # time data stored in SWOT h/w files
    SWOThwdir="C:/Users/coss.31/OneDrive - The Ohio State University/Documents/PYfun/valid/SWOTdirect/input/swot"
    onlyfiles = [f for f in listdir( SWOThwdir) if isfile(join( SWOThwdir, f))]
    IDstr=IDXnc[-25:-14]
    Tfile = [match for match in onlyfiles if IDstr in match]    
    ncf = Dataset(SWOThwdir+'/'+Tfile[0])
    T=np.ma.getdata(ncf.variables['nt'][:])
    SWt=[]
    for time in T:        
        d=datetime.strptime(str(time), "%Y%m%d")
        SWt.append(d.toordinal()+1)
    #pull q next    
    ncf = Dataset(IDXnc)
    #parse groups
    geobam=ncf.groups['geobam']
    hivdi=ncf.groups['hivdi']
    metroman=ncf.groups['metroman']
    momma=ncf.groups['momma']
    #put Qs into single SWq dict
    SWq={
        "geobam":np.ma.getdata(geobam.variables['q'][:]),
        "hivdi":np.ma.getdata(hivdi.variables['q'][:]),
        "metroman":np.ma.getdata(metroman.variables['q'][:]),
        "momma":np.ma.getdata(momma.variables['q'][:])
        }
    return SWt,SWq
     
    

    



 