# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 13:31:36 2021

@author: coss.31
@author: efried130
"""
# Standard imports
import datetime

# Third-party imports
import matplotlib.pyplot as plt  
import matplotlib.dates as Pdate
import numpy as np
from scipy.stats import pearsonr


def stats(St, Sq_, Vt, Vq, gid, IDstr, figdir):
    EMPTY = -9999
    MIN_OBS = 20  # Minimum observations required for metric calculation
    
    OUTconsensus = Sq_.get('consensus', EMPTY)
    
    ## trim and align
    Ro = []
    Gido = []
    SIGeo = []
    NSEo = []
    Rsqo = []
    KGEo = []
    RMSEo = []
    nRMSEo = []
    nBIASo = []
    no = []
    offkey = []
    STST = []

    for algo in Sq_:
        print(algo)

        Sq = Sq_[algo]
        if np.size(Sq) > 3:  # prevent failures on filled algo data
            STfilter = ~np.isnan(St)
            goodst = np.array(St)[STfilter]
            goodflpeq = Sq[STfilter]

            gqindex = []
            sqindex = []
            for time in goodst:
                tindex = np.where(Vt.astype(int) == int(time))[0]
                Stindex = np.where(goodst == int(time))[0]
                if np.size(tindex) > 0:
                    gqindex.append(tindex[0])
                    sqindex.append(Stindex[0])
                else:
                    gqindex.append(np.nan)
                    sqindex.append(np.nan)

            Ggqindex = ~np.isnan(gqindex)
            GQdx = np.array(gqindex)[Ggqindex].astype(int)
            
            Gsqindex = ~np.isnan(sqindex)
            GQsdx = np.array(sqindex)[Gsqindex].astype(int)
            
            Filtergq = np.array(Vq)[GQdx]
            Filterst = goodst[GQsdx]
            Filtersoneq = goodflpeq[GQsdx]

            # Filter to positive and upper ceiling values on both sides
            finalfilter = (Filtersoneq > 0) & (Filtergq > 0) & (Filtersoneq < 1e7) & (Filtergq < 1e7)

            Gq = Filtergq[finalfilter]
            St_p = Filterst[finalfilter]
            Sq = Filtersoneq[finalfilter]

            # Require minimum number of observations
            if len(Sq) >= MIN_OBS:

                mean_Gq = np.mean(Gq)

                # Pearson r
                r, _ = pearsonr(Sq, Gq)
                Ro.append(r)

                Gido.append(gid)
                offkey.append(algo)

                # sigmaE
                one_sigma = np.std(Gq - Sq) / mean_Gq if mean_Gq != 0 else EMPTY
                SIGeo.append(one_sigma)

                # NSE
                ss_res = np.sum((Gq - Sq) ** 2)
                ss_tot = np.sum((Gq - mean_Gq) ** 2)
                NSE_ = 1 - (ss_res / ss_tot) if ss_tot != 0 else EMPTY
                NSEo.append(NSE_)

                # Rsq (based on Pearson r)
                Rsq = r ** 2 if not np.isnan(r) else EMPTY
                Rsqo.append(Rsq)

                # KGE (based on Pearson r)
                if not np.isnan(r) and np.std(Gq) != 0 and mean_Gq != 0:
                    KGE_ = 1 - np.sqrt(
                        (r - 1) ** 2 +
                        ((np.std(Sq) / np.std(Gq)) - 1) ** 2 +
                        ((np.mean(Sq) / mean_Gq) - 1) ** 2
                    )
                else:
                    KGE_ = EMPTY
                KGEo.append(KGE_)

                # n
                n_ = len(Gq)
                no.append(n_)

                # RMSE
                RMSE = np.sqrt(np.mean((Sq - Gq) ** 2))
                RMSEo.append(RMSE)

                # nRMSE
                nRMSE = RMSE / mean_Gq if mean_Gq != 0 else EMPTY
                nRMSEo.append(nRMSE)

                # nBIAS
                BIAS = np.sum(Sq - Gq) / len(Gq)
                nBIAS = BIAS / mean_Gq if mean_Gq != 0 else EMPTY
                nBIASo.append(nBIAS)

                STST.append(St_p)

                # Time series plots
                if len(Sq) >= MIN_OBS:
                    strDates = []
                    for day in St_p:
                        date = datetime.date.fromordinal(int(day))
                        strDates.append(date)
                    fig = plt.figure()
                    ax = fig.add_subplot()
                    dates = Pdate.date2num(strDates)
                    ax.plot_date(dates, Gq, fmt='-', marker='o')
                    ax.plot_date(dates, Sq, fmt='-', marker='o')
                    ax.set_ylabel('Q (m^3/s)')
                    ax.legend(['Gage', algo])
                    fig.autofmt_xdate()
                    figname = f"{figdir}/{IDstr}_{algo}.jpg"
                    fig.savefig(figname)
                    plt.close(fig)

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

    St = np.array(St)
    St[np.isnan(St)] = EMPTY

    validout = {
        "algorithm": np.array([offkey]),
        "Gid": np.array(Gido),
        "pearsonr": Ro[:],
        "SIGe": SIGeo[:],
        "NSE": NSEo[:],
        "Rsq": Rsqo[:],
        "KGE": KGEo[:],
        "RMSE": RMSEo[:],
        "nRMSE": nRMSEo[:],
        "nBIAS": nBIASo[:],
        "n": no[:],
        "t": St,
        "consensus": OUTconsensus,
    }

    return validout