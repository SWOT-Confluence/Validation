# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 11:10:20 2021

@author: coss.31
"""

def write(Vout,IDstr,SWt):
   dump = "C:/Users/coss.31/OneDrive - The Ohio State University/Documents/PYfun/Valid/ValidDump/"
   EMPTY=-9999
   import netCDF4 as nc   
   from datetime import datetime
  
   now = datetime.now()
   now=now.strftime("%Y-%m-%d") 
  
   fn = dump+IDstr+'_validation.nc'
   ds = nc.Dataset(fn, 'w', format='NETCDF4') 
   ds.description= 'Statistics for reach_'+IDstr
   ds.history= 'File generated '+now
   if not Vout ==EMPTY:
       ds.hasvalidation='True'
       #dim
       Adim = ds.createDimension('Adim', len(Vout['algorithm']))  
       Tdim=ds.createDimension('Tdim', len(SWt))
       # create var
       algorithm = ds.createVariable('algorithim', 'str', ('Adim',))
      
       #
       NSE = ds.createVariable('NSE', 'f4', ('Adim',))  
       Rsq= ds.createVariable('Rsq', 'f4', ('Adim',))   
       KGE = ds.createVariable('KGE', 'f4', ('Adim',))  
       RMSE = ds.createVariable('RMSE', 'f4', ('Adim',))
       RMSE.units = 'm^3/s'   
       n = ds.createVariable('testn', 'f4', ('Adim',))
       #
       INdates= ds.createVariable('INdates', 'f4', ('Tdim',))
       INdates.units = 'days since Jan 1 Year 1'
       # fill var
       alglist=list(Vout['algorithm'])
       for alg in range(len(alglist)):
           algorithm[alg]= alglist[alg]    
       
        #
       NSE[:] = Vout['NSE']
       Rsq[:] = Vout['Rsq']
       KGE[:] = Vout['KGE']
       RMSE[:] = Vout['RMSE']
       n[:] = Vout['n']
       INdates[:]=SWt
   else:
       ds.hasvalidation='False'
       #dim
       Adim = ds.createDimension('Adim',1)  
       Tdim=ds.createDimension('Tdim',1)
       # create var
       algorithm = ds.createVariable('algorithim', 'S1', ('Adim',))
      
       #
       NSE = ds.createVariable('NSE', 'f4', ('Adim',))  
       Rsq= ds.createVariable('Rsq', 'f4', ('Adim',))   
       KGE = ds.createVariable('KGE', 'f4', ('Adim',))  
       RMSE = ds.createVariable('RMSE', 'f4', ('Adim',))
       RMSE.units = 'm^3/s'   
       n = ds.createVariable('testn', 'f4', ('Adim',))
       #
       INdates= ds.createVariable('INdates', 'f4', ('Tdim',))
       INdates.units = 'days since Jan 1 Year 1'
       # fill var  
       algorithm[:]= EMPTY
       #
       NSE[:] = EMPTY
       Rsq[:] = EMPTY
       KGE[:] = EMPTY
       RMSE[:] = EMPTY
       n[:] = EMPTY
       INdates[:]=EMPTY
  
   
   #close
   ds.close()