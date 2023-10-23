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


def stats(St,Sq,Vt,Vq,IDstr,figdir):
    
    # figdir="C:/Users/coss.31/OneDrive - The Ohio State University/Documents/PYfun/valid/figs"
    EMPTY=-9999
    St=np.array(St)
    # first define overlap times and filter timeseres before comparison
    #pull actual Q and t from uniform mostly empty block
    Vt_=Vt[np.logical_not(np.isnan(Vq))] 
    Vq_=Vq[np.logical_not(np.isnan(Vq))]
    #aditionally remove any -9999
    Vt_=Vt_[Vq_>0]
    Vq_=Vq_[Vq_>0]    
    olt,idV,ids=np.intersect1d(Vt_,St,return_indices=True)# these are still potential matches, need to filter again per SWOTq algorithim
    #t does not change between groups?
    Vt_=Vt_[idV]
    Vq_=Vq_[idV]
    # with nans out index overlap times
    NSEo=[]
    Rsqo=[]
    KGEo=[]
    RMSEo=[]
    nRMSEo=[]
    nBIASo=[]
    rRMSEo=[]
    no=[]
    offkey=[]
    if len(Sq)<10:
        for Grp in Sq:    
            Sq_=Sq[Grp]       
            Sq_=Sq_[ids]           
            offkey.append(Grp)
            ## strictly for testing I am converting -999999999999 to an array of posative value
            #Sq_=np.random.uniform(low=0.0, high=300.0,size=len(Sq_))
            if any(Sq_>EMPTY) and any(Vq_>EMPTY):
                 #SWOT empty removal
                St_=St[Sq_>0]
                Sq_=Sq_[Sq_>0]
                olt,idV,ids=np.intersect1d(Vt_,St_,return_indices=True)
                Vt_t=Vt_[idV]
                Vq_t=Vq_[idV]
                St_ = St_[ids]
                Sq_ = Sq_[ids]
                
                # NSE
                top=np.sum((Vq_t-Sq_)**2)
                bottom=np.sum(Vq_t-np.mean(Vq_t)**2)
                NSE=1-(top/bottom)
                NSEo.append(NSE)
                 #Rsq
                r=np.corrcoef( Sq_,Vq_t)
                r=r[0,1]
                Rsq=r**2
                Rsqo.append(Rsq)
                #KGE
                KGE=1-np.sqrt((r-1)**2 + ((np.std( Sq_)/np.std(Vq_t))-1)**2 +((np.mean(Sq_)/np.mean(Vq_t))-1)**2)  
                KGEo.append(KGE)
                 #n
                n=len(Vq_t)
                no.append(n)
                #RMSE
                RMSE=np.sqrt((np.sum( (Sq_ - Vq_t)**2))/n)
                RMSEo.append(RMSE)                
                #nRMSE
                NRMSE=RMSE/np.mean(Vq_t)
                nRMSEo.append(NRMSE)
                #nBIASo
                BIAS= np.sum(Sq_ - Vq_t)/len( Vq_t)
                nBIAS=BIAS/np.mean(Vq_t)
                nBIASo.append(nBIAS)
                #rRMSE
                rRMSEo.append(np.sqrt(NRMSE**2-nBIAS**2))
                
            else:
                    NSEo.append(EMPTY)
                    Rsqo.append(EMPTY)
                    KGEo.append(EMPTY)
                    no.append(EMPTY)
                    RMSEo.append(EMPTY)
                    nRMSEo.append(EMPTY)
                    nBIASo.append(EMPTY)
                    rRMSEo.append(EMPTY)                    
       
        validout={
            "algorithm": np.array([offkey]),
            "NSE":NSEo[:],
            "Rsq":Rsqo[:],
            "KGE":KGEo[:],
            "RMSE":RMSEo[:],
            "nRMSE":nRMSEo[:],
            "nBIAS":nBIASo[:],
            "rRMSE":rRMSEo[:],
            "n":no[:],
            "t":St
            }
    else:
        
        for Grp in Sq:
            Sq_=Sq[Grp]
            Sq_=Sq_[ids]           
            St__= St[ids]
            alg=Grp
            offkey.append(Grp)
            if any(Sq_>EMPTY) and any(Vq_>EMPTY):
                 #SWOT empty removal
                St_=St__[Sq_>0]
                Sq_=Sq_[Sq_>0]
                olt,idV,ids=np.intersect1d(Vt_,St_,return_indices=True)
                Vt_t=Vt_[idV]
                Vq_t=Vq_[idV]
                St_ = St_[ids]
                Sq_ = Sq_[ids]
                #plot and save
                strDates = []               
                for day in Vt_t:
                    date= datetime.date.fromordinal(day)
                    strDates.append(date)
                fig = plt.figure()                
                ax = fig.add_subplot()
                dates = Pdate.date2num(strDates)
                ax.plot_date(dates, Vq_t, fmt='-',marker='o')
                ax.plot_date(dates, Sq_, fmt='-',marker='o')               
               
                ax.set_ylabel('Q (m^3/s)')
                ax.legend(['Gage',alg])
                fig.autofmt_xdate()
                # figname=figdir+'/'+IDstr+alg+'.jpg'
                figname=f"{figdir}/{IDstr}_{alg}.jpg"
                fig.savefig(figname)
                
                # NSE
                top=np.sum((Vq_t-Sq_)**2)
                bottom=np.sum(Vq_t-np.mean(Vq_t)**2)
                NSE=1-(top/bottom)
                NSEo.append(NSE)
                 #Rsq
                r=np.corrcoef( Sq_,Vq_t)
                r=r[0,1]
                Rsq=r**2
                Rsqo.append(Rsq)
                #KGE
                KGE=1-np.sqrt((r-1)**2 + ((np.std( Sq_)/np.std(Vq_t))-1)**2 +((np.mean(Sq_)/np.mean(Vq_t))-1)**2)  
                KGEo.append(KGE)
                 #n
                n=len(Vq_t)
                no.append(n)
                #RMSE
                RMSE=np.sqrt((np.sum( (Sq_ - Vq_t)**2))/n)
                RMSEo.append(RMSE)                
                #nRMSE
                NRMSE=RMSE/np.mean(Vq_t)
                nRMSEo.append(NRMSE)
                #nBIASo
                BIAS= np.sum(Sq_ - Vq_t)/len( Vq_t)
                nBIAS=BIAS/np.mean(Vq_t)
                nBIASo.append(nBIAS)
                #rRMSE
                rRMSEo.append(np.sqrt(NRMSE**2-nBIAS**2))
            else:
                    NSEo.append(EMPTY)
                    Rsqo.append(EMPTY)
                    KGEo.append(EMPTY)
                    no.append(EMPTY)
                    RMSEo.append(EMPTY)
                    nRMSEo.append(EMPTY)
                    nBIASo.append(EMPTY)
                    rRMSEo.append(EMPTY) 
        validout={
            "algorithm": np.array([offkey]),
            "NSE":NSEo[:],
            "Rsq":Rsqo[:],
            "KGE":KGEo[:],
            "RMSE":RMSEo[:],
            "nRMSE":nRMSEo[:],
            "nBIAS":nBIASo[:],
            "rRMSE":rRMSEo[:],
            "n":no[:],
            "t":St
            }
   
    return  validout
    
