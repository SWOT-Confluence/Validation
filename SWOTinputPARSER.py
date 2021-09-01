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
        #Q
        USGSq = np.ma.getdata(ncf.variables['USGS_Q'][:])    
        #Qt
        USGSt = np.ma.getdata(ncf.variables['USGS_Qt'][:]) 
        return  USGS_reach, USGSq,  USGSt
    elif 'GRDC' in IDXnc:
         # need to pull Q, T, and reach identifiers
        #reach
        GRDC_reach=np.ma.getdata(ncf.variables['Reach_ID'][:])   
        #Q
        GRDCq = np.ma.getdata(ncf.variables['GRDC_Q'][:])    
        #Qt
        GRDCt = np.ma.getdata(ncf.variables['GRDC_Qt'][:]) 
            
        return GRDC_reach, GRDCq,  GRDCt

    



 