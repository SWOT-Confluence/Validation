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
    # need to 
    dataUSGS=np.ma.getdata(ncf.variables['STAID'][:])
    dataUSGS=np.char.decode(dataUSGS)
    #
    reachID = np.ma.getdata(ncf.variables['reach_id'][:])
    reachID=np.char.decode(reachID)
    USt={}
    Rt={}
    for i in range(len(reachID)):
        US=','.join(dataUSGS[i,:])
        USt[i]=US.replace(',','')
        #
        R=','.join( reachID[i,:])
        Rt[i]=R.replace(',','')
        
    dataUSGS = USt
    reachID = Rt
    return dataUSGS, reachID
def flag(In): 
    In = In.replace(np.nan,'*', regex=True)
   


# e   Value has been edited or estimated by USGS personnel and is write protected
# &     Value was computed from affected unit values
# E     Value was computed from estimated unit values.
# A     Approved for publication -- Processing and review completed.
# P     Provisional data subject to revision.
# <     The value is known to be less than reported value and is write protected.
# >     The value is known to be greater than reported value and is write protected.
# 1     Value is write protected without any remark code to be printed
# 2     Remark is write protected without any remark code to be printed
#       No remark (blank)
    M={}
    for i in range(len(In)):
        if 'Ice' in In[i] or '*' in In[i]:
            M[i]=False
        else:
            M[i]=True
                
                
    Mask=pd.array(list(M.values()),dtype="boolean")
    return Mask
 