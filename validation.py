# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 13:31:36 2021

@author: coss.31
"""

def stats(St,Sq,Vt,Vq):
    import numpy as np
    EMPTY=-9999
    # first define overlap times and filter timeseres before comparison
    #pull actual Q and t from uniform mostly empty block
    Vt_=Vt[np.logical_not(np.isnan(Vq))] 
    Vq_=Vq[np.logical_not(np.isnan(Vq))]
    olt,idV,ids=np.intersect1d(Vt_,St,return_indices=True)
    #t does not change between groups?
    Vt_=Vt_[idV]
    Vq_=Vq_[idV]
    # with nans out index overlap times
    NSEo=[]
    Rsqo=[]
    KGEo=[]
    RMSEo=[]
    no=[]
    for Grp in Sq:    
        Sq_=Sq[Grp]       
        Sq_=Sq_[ids]
        ## strictly for testing I am converting -999999999999 to an array of posative value
        #Sq_=np.random.uniform(low=0.0, high=300.0,size=len(Sq_))
            
        if any(Sq_>EMPTY) and any(Vq_>EMPTY):     
           
            # NSE
            top=np.sum((Sq_-Vq_)**2)
            bottom=np.sum(Vq_-np.mean(Vq_))**2
            NSE=1-(top/bottom)
            NSEo.append(NSE)
             #Rsq
            r=np.corrcoef( Sq_,Vq_)
            r=r[0,1]
            Rsq=r**2
            Rsqo.append(Rsq)
            #KGE
            KGE=1-np.sqrt((r-1)**2 + ((np.std( Sq_)/np.std(Vq_))-1)**2 +((np.mean(Sq_)/np.mean(Vq_))-1)**2)  
            KGEo.append(KGE)
             #n
            n=len(Vq_)
            no.append(n)
            #RMSE
            RMSE=np.sqrt((np.sum( Sq_ - Vq_)**2)/n)
            RMSEo.append(RMSE)
        else:
                NSEo.append(EMPTY)
                Rsqo.append(EMPTY)
                KGEo.append(EMPTY)
                no.append(EMPTY)
                RMSEo.append(EMPTY)
    validout={
        "algorithm": ['geobam','hivdi','metroman','momma'],
        "NSE":NSEo[:],
        "Rsq":Rsqo[:],
        "KGE":KGEo[:],
        "RMSE":RMSEo[:],
        "n":no[:]      
        }
   
    return  validout
    