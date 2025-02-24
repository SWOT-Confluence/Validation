# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 13:31:36 2021

@author: coss.31
"""
# Standard imports
import datetime

# Third-party imports
import matplotlib.pyplot as plt  
import matplotlib.dates as Pdate
import numpy as np
from scipy.stats import spearmanr


def stats(St,Sq_,Vt,Vq,gid,IDstr,figdir):
    EMPTY=-9999    
    ## trim and allign
    Ro=[]
    Gido=[]
    SIGeo=[]
    NSEo=[]
    Rsqo=[]
    KGEo=[]
    RMSEo=[]
    nRMSEo=[]
    nBIASo=[]   
    no=[]
    offkey=[]
    STST=[]
    for algo in Sq_:
        print(algo)
        if algo =='consensus':
            OUTconsensus=Sq_[algo]
        else:
           OUTconsensus=EMPTY
       
        Sq=Sq_[algo]
        if np.size(Sq)>3:#this will prevent failures on filled algo data         
            STfilter=~np.isnan(St)
            goodst=np.array(St)[STfilter]
            
            goodflpeq=Sq[STfilter]  
            gqindex=[]
            sqindex=[]
            for time in goodst:
                    tindex=np.where(Vt.astype(int)==int(time))[0]
                    Stindex=np.where(goodst==int(time))[0]
                    if np.size(tindex)>0:
                        gqindex.append(tindex[0])
                        sqindex.append(Stindex[0])
                    else:
                        gqindex.append(np.nan)
                        sqindex.append(np.nan)

            Ggqindex=~np.isnan(gqindex)
            GQdx=np.array(gqindex)[Ggqindex]
            GQdx=GQdx.astype(int)

            Gsqindex=~np.isnan(sqindex)
            GQsdx=np.array(sqindex)[Gsqindex]
            GQsdx=GQsdx.astype(int)

            Filtergq=np.array(Vq)[GQdx]
            Filterst=goodst[GQsdx]
            Filtersoneq=goodflpeq[GQsdx]    
            flpefinalfilter=Filtersoneq>0
            gfinalfilter=Filtergq>0
            finalfilter= flpefinalfilter&gfinalfilter



            Gq=Filtergq[finalfilter]
            St_p=Filterst[finalfilter]    
            Sq=Filtersoneq[finalfilter]
            # Do the Stats
            if len(Sq)>3: #we need three values to not break things with some of these stats functions
               
                res= spearmanr(Sq, Gq)
                r=res.statistic
                Ro.append(r)       
                Gido.append(gid)     
                offkey.append(algo)            
                #sigmaE
                SIGeo.append(np.std((Gq-Sq))/np.mean(Gq))
                # NSE
                top=np.sum((Gq-Sq)**2)
                bottom=np.sum((Gq-np.mean(Gq))**2)
                NSE_=1-(top/bottom)
                NSEo.append(NSE_)
                #Rsq        
                Rsq=r**2
                Rsqo.append(Rsq)
                #KGE
                KGE_=1-np.sqrt((r-1)**2 + ((np.std( Sq)/np.std(Gq))-1)**2 +((np.mean(Sq)/np.mean(Gq))-1)**2)
                KGEo.append(KGE_)
                #n
                n_=len(Gq)
                no.append(n_)
                #RMSE
                RMSE=np.sqrt((np.sum( (Sq - Gq)**2))/n_)
                RMSEo.append(RMSE)
                #nRMSE
                NRMSE=RMSE/np.mean(Gq)
                nRMSEo.append(NRMSE)
                #nBIASo
                BIAS= np.sum(Sq - Gq)/len(Gq)
                nBIAS=BIAS/np.mean(Gq)
                nBIASo.append(nBIAS)
                STST.append(St_p)
                if len(Sq)>10: #We still want all our ts plots
                    #plot and save
                    strDates = []
                 
                    for day in St_p:
                        date= datetime.date.fromordinal(int(day))
                        strDates.append(date)
                    fig = plt.figure()                
                    ax = fig.add_subplot()
                    dates = Pdate.date2num(strDates)
                    ax.plot_date(dates, Gq, fmt='-',marker='o')
                    ax.plot_date(dates, Sq, fmt='-',marker='o')               

                    ax.set_ylabel('Q (m^3/s)')
                    ax.legend(['Gage',algo])
                    fig.autofmt_xdate()                
                    figname=f"{figdir}/{IDstr}_{algo}.jpg"
                    fig.savefig(figname)      
       
            else:
                Ro.append(EMPTY)
                Gido.append(gid)
                SIGeo.append(EMPTY)
                NSEo.append(EMPTY)
                Rsqo.append(EMPTY)
                KGEo.append(EMPTY)
                RMSEo.append(EMPTY) 
                nRMSEo.append(EMPTY)
                nBIASo.append(EMPTY)      
                no.append(len(Sq))
                offkey.append(algo)
                STST.append(EMPTY)

        else:
            Ro.append(EMPTY)
            Gido.append(gid)
            SIGeo.append(EMPTY)
            NSEo.append(EMPTY)
            Rsqo.append(EMPTY)
            KGEo.append(EMPTY)
            RMSEo.append(EMPTY) 
            nRMSEo.append(EMPTY)
            nBIASo.append(EMPTY)      
            no.append(EMPTY)
            offkey.append(algo)
            STST.append(EMPTY)
    
    St=np.array(St)
    St[np.isnan(St)]=EMPTY   
    validout={
        "algorithm": np.array([offkey]),
        "Gid": np.array(Gido),
        "Spearmanr":Ro[:],
        "SIGe":SIGeo[:],
        "NSE":NSEo[:],
        "Rsq":Rsqo[:],
        "KGE":KGEo[:],
        "RMSE":RMSEo[:],
        "nRMSE":nRMSEo[:],
        "nBIAS":nBIASo[:],            
        "n":no[:],
        "t": St,
        "consensus":OUTconsensus
        }
   
    return  validout
    
