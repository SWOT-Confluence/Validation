"""Module to run validation operations and output stats.

Runs on a reach of data and requires JSON data for reach retrieved by
AWS Batch index.

Class
-----
ValidationConfluence: Stores data and executes validation operations.

Constants
---------
INPUT_DIR: Path
    path to input directory
OFFLINE_DIR: Path
    path to offline directory
OUTPUT_DIR: Path
    path to output directory

Functions
---------
get_reach_data(input_json)
    return dictionary of reach data
run_validation()
    orchestrate validation operations
"""

# Standard imports
import datetime
import json
import os
from pathlib import Path
import sys
import warnings
import matplotlib.pyplot as plt
import seaborn as sb

# Local imports
from val.validation import stats

# Third-party imports
from netCDF4 import Dataset, stringtochar,chartostring
import numpy as np

# Constants
INPUT = Path("/mnt/data/input")
FLPE=Path("/mnt/data/flpe")
MOI=Path("/mnt/data/moi")
OFFLINE = Path("/mnt/data/offline")
OUTPUT = Path("/mnt/data/output")
# INPUT = Path('/Users/mtd/Analysis/SWOT/Discharge/Confluence/verify/validation/input')
# OFFLINE = Path('/Users/mtd/Analysis/SWOT/Discharge/Confluence/verify/validation/offline')
# OUTPUT = Path('/Users/mtd/Analysis/SWOT/Discharge/Confluence/verify/validation/output')

class ValidationConfluence:
    """Class that runs validation operations for Confluence workflow.
    
    Attributes
    ----------
    gage_data: dict
        dictionary of gage reach identifiers, q, and qt
    input_dir: Path
        path to input directory
    INT_FILL: int
        integer fill value used in NetCDF files
    NUM_ALGOS: int
        number of algorithms to store data for
    offline_data: dict
        dictionary of offline discharge values and time
    output_dir: Path
        path to output directory
    reach_id: int
        unique reach identifier
    
    Methods
    -------
    read_gage_data()
        read gage data from SoS file
    get_gage_q(sos, gage_type)
        return discharge and discharge time from gage_type
    is_offline_valid(offline_data)
        check if offline data is only comprised of NaN values
    read_offline_data(reach_id)
        reads data from offline module and stores in flpe_data dictionary
    is_flpe_valid(flpe_data)
        check if flpe data is only comprised of NaN values
    read_flpe_data(reach_id)
        reads data from flpe module and stores in flpe_data dictionary
    is_moi_valid(moi_data)
        check if moi data is only comprised of NaN values
    read_moi_data(reach_id)
        reads data from moi module and stores in flpe_data dictionary
    read_time_data()
        read time of observations from SWOT files
    validate()
        run validation operations on gage data and FLPE data; write stats
    write(stats, time, reach_id, gage_type)
        write stats to NetCDF file
    """

    INT_FILL = -999
    NUM_ALGOS = 14

    def __init__(self, reach_data,flpe_dir,moi_dir, offline_dir, input_dir, output_dir, run_type):

        """
        Parameters
        ----------
        reach_data: dict
            dictionary of reach identifier and associated file names
        offline_dir: Path
            path to offline data directory
        input_dir: Path 
            path to input directory
        output_dir: Path
            path to output directory
        run_type: str
            string indicating if we are doing a constrained or unconstrained run
        """
        
        self.input_dir = input_dir
        self.run_type = run_type
        self.reach_id = reach_data["reach_id"]
        print('Processing', self.reach_id)
        self.gage_data = self.read_gage_data(input_dir / "sos" / reach_data["sos"])
        self.offline_data = self.read_offline_data(offline_dir)
        self.flpe_data = self.read_flpe_data(flpe_dir)
        self.moi_data = self.read_moi_data(moi_dir)
        self.output_dir = output_dir


    def read_gage_data(self, sos_file):
        """Read gage data from SoS file and stores in gage data dictionary."""

        sos = Dataset(sos_file, 'r')
        gage_data = {}
        # could be optimized using the global gage_agency varaible
        for gage_agency in sos.gauge_agency.split(';'):
            gage_data = self.get_gage_q(sos, gage_agency)
            if gage_data != {}:
                print('found gage')
                break

        # groups = list(sos.groups.keys())
        # gage_data = {}
        # print('finding groups')
        
        # if "USGS" in groups:
        #     print('in usgs')
        #     gage_data = self.get_gage_q(sos, "USGS")
        #     if gage_data == {}:
        #         gage_data = self.get_gage_q(sos, "WSC")
        # elif "grdc" in groups:
        #     gage_data = self.get_gage_q(sos, "grdc")
        # elif "Hidroweb" in groups:
        #     gage_data = self.get_gage_q(sos, "Hidroweb")
        # elif "ABOM" in groups:
        #     gage_data = self.get_gage_q(sos, "ABOM")
        # elif "MLIT" in groups:
        #     gage_data = self.get_gage_q(sos, "MLIT")
        # elif "DEFRA" in groups:
        #     gage_data = self.get_gage_q(sos, "DEFRA")
        #     if gage_data == {}:
        #         gage_data = self.get_gage_q(sos, "EAU")
            
        sos.close()
        return gage_data
    
    def get_gage_q(self, sos, gage_type):
        """Return discharge and discharge time from gage_type.
        
        gage_type values should be either 'usgs' or 'grdc'.
        
        Parameters
        ----------
        sos: NetCDF dataset
            SOS NetCDF dataset reference
        gage_type: str
            indicates type of gage to search for
        
        Returns
        -------
        dictionary of discharge and discharge time
        """
        
        gage = sos[gage_type]
        rids = gage[f"{gage_type}_reach_id"][:].filled(np.nan)
        index = np.where(self.reach_id == rids)
        print('here is index we are workign with ', index)
        # if its more than one, we take it down to a scalar
        if len(index[0])>1:
             warnings.warn('multiple gages for this reach. Selecting closest meanQ to model')
             #pull model q for this reach
             modelindex=np.where(self.reach_id == sos['reaches']['reach_id'][:].filled(np.nan))
             model_q=sos['model']['mean_q'][modelindex][:].filled(np.nan)
             gmq=[]
             glt=[]
             for Gindex in index:
                 #pull mean q and timeseries lenghts
                 gmq.append(gage[f"{gage_type}_mean_q"][Gindex][:].filled(np.nan))
                 t = gage[f"{gage_type}_qt"][Gindex][:].filled(self.INT_FILL).astype(int)
                 glt.append(len(t[t>0]))
                 
             if np.isnan(model_q):
                 #when model is nan, choose longest timeseries
                 index=np.array(index[np.argmax(np.array(glt))])
                 if np.size(index)>1:
                     warnings.warn('model was nan and times are same length')
                     index=index[0]
                 
             else:
                 #othewise closest mean
                 index=np.array(index[np.argmin(np.abs(np.array(glt)-model_q))])
                 if np.size(index)>1:
                     warnings.warn('identical mean q values')
                     index=index[0]
        elif len(index[0]) == 1:
            index = index[0][0]
        

            

        gage_data = {}
        if np.isscalar(index):
            if self.run_type == "constrained":
                # if constraind check and see if the gage selected at this index is a 0
                if gage["CAL"][:][index] == 1:
                    warnings.warn('gauge found was calibration.. This is a constrained run, so it will not be used for validation')
                    return gage_data
            gage_data["type"] = gage_type
            gage_data["q"] = gage[f"{gage_type}_q"][index][:].filled(np.nan)
            gage_data["qt"] = gage[f"{gage_type}_qt"][index][:].filled(self.INT_FILL).astype(int)
            gage_data["gid"] = chartostring(gage[f"{gage_type}_id"][index][:].filled(np.nan))
            
        return gage_data
    
    def read_moi_data(self,moi_dir):
        """Reads data from moi module and returns dictionary.
        
        Parameters
        ----------
        moi_dir: Path
            path to moi data directory
        
        Returns
        -------
        dictionary of algorithm moi results
        """
       

        moi_file = f"{moi_dir}/{self.reach_id}_integrator.nc"
        moi = Dataset(moi_file, 'r')
        moi_data = {}
        moi_data["metroman"] =  moi["metroman/q"][:].filled(np.nan)
        moi_data["neobam"] =  moi["neobam/q"][:].filled(np.nan)
        moi_data["hivdi"] =  moi["hivdi/q"][:].filled(np.nan)
        moi_data["momma"] =  moi["momma/q"][:].filled(np.nan)
        moi_data["sad"] = moi["sad/q"][:].filled(np.nan)
        moi_data["sic4dvar"] = moi["sic4dvar/q"][:].filled(np.nan)      
        moi.close()      
        
        #create pre-offline consensus
        ALLQ=np.full((len(moi_data.keys()), len(moi_data["metroman"])), np.nan)
        for row in range(len(moi_data.keys())):
            ALGv=moi_data[list(moi_data.keys())[row]]        
            ALGv[ALGv<0]=np.nan
            ALLQ[row,:]=ALGv

        consensus=np.nanmedian(ALLQ,axis=0)
        moi_data["consensus"]=consensus
        
        
        if self.is_moi_valid(moi_data):
            return moi_data
        else: 
            return {}
    
    def is_moi_valid(self, moi_data):
        """Check if moi data is only comprised of NaN values.
        
        Returns
        -------
        False if all NaN values are present otherwise True
        """
        
        invalid = 0
        for v in moi_data.values():
            if np.count_nonzero(~np.isnan(v)) == 0: invalid += 1
        if invalid == self.NUM_ALGOS:
            print('MOI IS NOT VALID')
            return False
        else:
            return True
        
    def read_flpe_data(self,flpe_dir):
        """Reads data from flpe module and returns dictionary.
        
        Parameters
        ----------
        flpe_dir: Path
            path to flpe data directory
        
        Returns
        -------
        dictionary of algorithm moi results
        """
        convention_dict = {
            "metroman":"allq",
            "neobam":"q/q",
            "hivdi":"reach/Q",
            "momma":"Q",
            "sad":"Qa",  
            "sic4dvar":"Q_da",        
            
            
        }

        flpe_file_metroman = f"{flpe_dir}/{'metroman'}/{self.reach_id}_metroman.nc"
        flpe_file_neobam = f"{flpe_dir}/{'geobam'}/{self.reach_id}_geobam.nc"
        flpe_file_hivdi = f"{flpe_dir}/{'hivdi'}/{self.reach_id}_hivdi.nc"
        flpe_file_momma = f"{flpe_dir}/{'momma'}/{self.reach_id}_momma.nc"
        flpe_file_sad = f"{flpe_dir}/{'sad'}/{self.reach_id}_sad.nc"
        flpe_file_sic4dvar = f"{flpe_dir}/{'sic4dvar'}/{self.reach_id}_sic4dvar.nc"
        try:
            flpe_mm = Dataset(flpe_file_metroman, 'r')
        except:
            flpe_mm=-9999
        try:    
            flpe_nb = Dataset(flpe_file_neobam, 'r')
        except:
            flpe_nb=-9999
        try:
            flpe_hi = Dataset(flpe_file_hivdi, 'r')
        except:
            flpe_hi=-9999
        try:
            flpe_mo = Dataset(flpe_file_momma, 'r')
        except:
            flpe_mo=-9999
        try:
            flpe_sa = Dataset(flpe_file_sad, 'r')
        except:
            flpe_sa=-9999
        try:
            flpe_si = Dataset(flpe_file_sic4dvar, 'r')
        except:
            flpe_si=-9999
      
        
        flpe_data = {}
        
        if flpe_mm==-9999:
            flpe_data["metroman"]=-9999
        else:
            flpe_data["metroman"] =  flpe_mm[convention_dict["metroman"]][:].filled(np.nan)
            conlen=len(flpe_data["metroman"])
            flpe_mm.close()
        if flpe_nb==-9999:
            flpe_data["neobam"] =-9999
        else:    
            flpe_data["neobam"] =  flpe_nb[convention_dict["neobam"]][:].filled(np.nan)
            conlen=len(flpe_data["neobam"])
            flpe_nb.close()
        if flpe_hi ==-9999:
            flpe_data["hivdi"] =-9999
        else:
            flpe_data["hivdi"] =  flpe_hi[convention_dict["hivdi"]][:].filled(np.nan)
            conlen=len(flpe_data["hivdi"])
            flpe_hi.close()
        if  flpe_mo==-9999:
            flpe_data["momma"]=-9999
        else:            
            flpe_data["momma"] =  flpe_mo[convention_dict["momma"]][:].filled(np.nan)
            conlen=len(flpe_data["momma"])
            flpe_mo.close()
        if flpe_sa==-9999:
            flpe_data["sad"]=-9999
        else:
            flpe_data["sad"] = flpe_sa[convention_dict["sad"]][:].filled(np.nan)
            conlen=len(flpe_data["sad"])
            flpe_sa.close()
        if flpe_si==-9999:
            flpe_data["sic4dvar"]=-9999
        else:
            flpe_data["sic4dvar"] = flpe_si[convention_dict["sic4dvar"]][:].filled(np.nan)
            conlen=len(flpe_data["sic4dvar"])
            flpe_si.close()

       
        #create pre-offline consensus
        ALLQ=np.full((len(flpe_data.keys()),  conlen), np.nan)
        for row in range(len(flpe_data.keys())):
            ALGv=flpe_data[list(flpe_data.keys())[row]]
            if np.size(ALGv)==conlen:
                ALGv[ALGv<0]=np.nan
                ALLQ[row,:]=ALGv
            

        consensus=np.nanmedian(ALLQ,axis=0)
        flpe_data["consensus"]=consensus
        
        
        if self.is_flpe_valid(flpe_data):
            return flpe_data
        else: 
            return {}
    
    def is_flpe_valid(self, flpe_data):
        """Check if moi data is only comprised of NaN values.
        
        Returns
        -------
        False if all NaN values are present otherwise True
        """
        
        invalid = 0
        for v in flpe_data.values():
            if np.count_nonzero(~np.isnan(v)) == 0: invalid += 1
        if invalid == self.NUM_ALGOS:
            print('flpe IS NOT VALID')
            return False
        else:
            return True                

    def read_offline_data(self, offline_dir):
        """Reads data from offline module and returns dictionary.
        
        Parameters
        ----------
        offline_dir: Path
            path to offline data directory
        
        Returns
        -------
        dictionary of algorithm offline results
        """
        convention_dict = {
            "metro_q_c":"dschg_gm",
            "bam_q_c":"dschg_gb",
            "hivdi_q_c":"dschg_gh",
            "momma_q_c":"dschg_go",
            "sads_q_c":"dschg_gs",
            "consensus_q_c":"dschg_gc",
            "sic4dvar_q_c":"dschg_gi",
            "metro_q_uc":"dschg_m",
            "bam_q_uc":"dschg_b",
            "hivdi_q_uc":"dschg_h",
            "momma_q_uc":"dschg_o",
            "sads_q_uc":"dschg_s",
            "sic4dvar_q_uc":"dschg_i",
            "consensus_q_uc":"dschg_c",
            "d_x_area":"d_x_area",
            "d_x_area_u":"d_x_area_u",
        }

        offline_file = f"{offline_dir}/{self.reach_id}_offline.nc"
        off = Dataset(offline_file, 'r')
        offline_data = {}
        offline_data[convention_dict["bam_q_c"]] = off[convention_dict["bam_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["hivdi_q_c"]] = off[convention_dict["hivdi_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["metro_q_c"]] = off[convention_dict["metro_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["momma_q_c"]] = off[convention_dict["momma_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["sads_q_c"]] = off[convention_dict["sads_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["sic4dvar_q_c"]] = off[convention_dict["sic4dvar_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["sic4dvar_q_uc"]] = off[convention_dict["sic4dvar_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["bam_q_uc"]] = off[convention_dict["bam_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["hivdi_q_uc"]] = off[convention_dict["hivdi_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["metro_q_uc"]] = off[convention_dict["metro_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["momma_q_uc"]] = off[convention_dict["momma_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["sads_q_uc"]] = off[convention_dict["sads_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["consensus_q_c"]] = off[convention_dict["consensus_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["consensus_q_uc"]] = off[convention_dict["consensus_q_uc"]][:].filled(np.nan)
        off.close()
        
        if self.is_offline_valid(offline_data):
            return offline_data
        else: 
            return {}
    
    def is_offline_valid(self, offline_data):
        """Check if offline data is only comprised of NaN values.
        
        Returns
        -------
        False if all NaN values are present otherwise True
        """
        
        invalid = 0
        for v in offline_data.values():
            if np.count_nonzero(~np.isnan(v)) == 0: invalid += 1
        if invalid == self.NUM_ALGOS:
            print('OFFLINE IS NOT VALID')
            return False
        else:
            return True

    def read_time_data(self):
        """Read time of observations from SWOT files.

        Parameters
        ----------
        reach_id: int
            unique reach identifier
        
        Returns
        -------
        list of ordinal times
        """

        swot = Dataset(self.input_dir / "swot" / f"{self.reach_id}_SWOT.nc", 'r')
        time = swot["reach"]["time"][:].filled(np.nan)
        swot.close()
        epoch = datetime.datetime(2000,1,1,0,0,0)
        ordinal_times = []

        for t in time:
            try:
                ordinal_times.append((epoch + datetime.timedelta(seconds=t)).toordinal())
            except:
                ordinal_times.append(np.nan)
                #print(time)
                warnings.warn('problem with time conversion to ordinal, most likely nan value')
        # return [ (epoch + datetime.timedelta(seconds=t)).toordinal() for t in time ]   # Check if this format works
        return ordinal_times

    def validate(self):
        """Run validation operations on gage data and FLPE data; write stats."""
     
        # SWOT time 
        time = self.read_time_data()
        algo_dim = int(self.NUM_ALGOS/2)
        # Data fill values
        data_flpe = {
            "algorithm": np.full( algo_dim, fill_value=""),
            "Gid": np.full(algo_dim, fill_value=""),
            "Spearmanr": np.full(algo_dim, fill_value=-9999),
            "SIGe": np.full(algo_dim, fill_value=-9999),
            "NSE": np.full(algo_dim, fill_value=-9999),           
            "Rsq": np.full(algo_dim, fill_value=-9999),
            "KGE": np.full(algo_dim, fill_value=-9999),           
            "RMSE": np.full(algo_dim, fill_value=-9999),           
            "n": np.full(algo_dim, fill_value=-9999),           
            "nRMSE":np.full(algo_dim, fill_value=-9999),           
            "nBIAS":np.full(algo_dim, fill_value=-9999),
            "t":np.full((self.NUM_ALGOS), fill_value=-9999),
            
           
        }

        no_flpe = False
        # Check if there is data to validate
        if self.gage_data:
            if self.flpe_data:
                #### Check should go here for all nan gauge data ---------------------------------
                data_flpe = stats(time, self.flpe_data, self.gage_data["qt"], 
                            self.gage_data["q"],self.gage_data["gid"], str(self.reach_id), 
                            self.output_dir / "figs")
            else:
                warnings.warn('No flpe data found...')
                no_flpe = True
        else:
            warnings.warn('No gauge found for reach...')

        
       
        data_moi = {
            "algorithm": np.full( algo_dim, fill_value=""),
            "Gid": np.full( algo_dim, fill_value=""),
            "Spearmanr": np.full( algo_dim, fill_value=-9999),
            "SIGe": np.full( algo_dim, fill_value=-9999),
            "NSE": np.full( algo_dim, fill_value=-9999),           
            "Rsq": np.full( algo_dim, fill_value=-9999),
            "KGE": np.full( algo_dim, fill_value=-9999),           
            "RMSE": np.full( algo_dim, fill_value=-9999),           
            "n": np.full( algo_dim, fill_value=-9999),           
            "nRMSE":np.full( algo_dim, fill_value=-9999),           
            "nBIAS":np.full( algo_dim, fill_value=-9999),
             "t":np.full((self.NUM_ALGOS), fill_value=-9999),
        }

        no_moi = False
        # Check if there is data to validate
        if self.gage_data:
            if self.moi_data:
                #### Check should go here for all nan gauge data ---------------------------------
                data_moi = stats(time, self.moi_data, self.gage_data["qt"], 
                            self.gage_data["q"],self.gage_data["gid"], str(self.reach_id), 
                            self.output_dir / "figs")
            else:
                warnings.warn('No moi data found...')
                no_moi = True
        else:
            warnings.warn('No gauge found for reach...')

        

        data_O = {
            "algorithm": np.full((self.NUM_ALGOS), fill_value=""),
            "Gid": np.full((self.NUM_ALGOS), fill_value=""),
            "Spearmanr": np.full((self.NUM_ALGOS), fill_value=-9999),
            "SIGe": np.full((self.NUM_ALGOS), fill_value=-9999),
            "NSE": np.full((self.NUM_ALGOS), fill_value=-9999),           
            "Rsq": np.full((self.NUM_ALGOS), fill_value=-9999),
            "KGE": np.full((self.NUM_ALGOS), fill_value=-9999),           
            "RMSE": np.full((self.NUM_ALGOS), fill_value=-9999),           
            "n": np.full((self.NUM_ALGOS), fill_value=-9999),           
            "nRMSE":np.full((self.NUM_ALGOS), fill_value=-9999),           
            "nBIAS":np.full((self.NUM_ALGOS), fill_value=-9999),
            "t":np.full((self.NUM_ALGOS), fill_value=-9999),
            
        }
        
        no_offline = False
        # Check if there is data to validate
        if self.gage_data:
            if self.offline_data:
                #### Check should go here for all nan gauge data ---------------------------------
                data_O = stats(time, self.offline_data, self.gage_data["qt"], 
                            self.gage_data["q"],  self.gage_data["gid"], str(self.reach_id), 
                            self.output_dir / "figs")
            else:
                warnings.warn('No offline data found...')
                no_offline = True
        else:
            warnings.warn('No gauge found for reach...')
            
        # Write out valid or invalid data
        gage_type = "No data" if not self.gage_data else self.gage_data["type"]       
        ALLnone=np.all([no_flpe,no_moi,no_offline])
       
        if (gage_type != "No data") and (ALLnone != True):
            self.write(data_flpe,data_moi,data_O, self.reach_id, gage_type,[no_flpe,no_moi,no_offline])

    def write(self, stats_flpe,stats_moi,stats_O, reach_id, gage_type,GO):
        FLPEno=GO[0]
        MOIno=GO[1]
        OFFno=GO[2]
        #print(stats_flpe)
        """Write stats to NetCDF file.
        
        Parameters
        ----------
        stats_flpe: dict
            dictionary of flpe stats for each algorithm
        stats_moi: dict
            dictionary of moi stats for each algorithm
        stats_O: dict
            dictionary of offline stats for each algorithm
        reach_id: int
            reach identifier for stats
        gage_type: str
            type of gage data used for validation
        """

        fill = -999999999999
        empty = -9999

        out = Dataset(self.output_dir / "stats" / f"{reach_id}_validation.nc", 'w')
        out.reach_id = reach_id
        out.description = f"Statistics for reach: {reach_id}"
        out.history = datetime.datetime.now().strftime('%d-%b-%Y %H:%M:%S')
        out.has_validation_flpe = 0 if np.where(stats_flpe["algorithm"] == "")[0].size == self.NUM_ALGOS/2 else 1
        out.has_validation_moi = 0 if np.where(stats_moi["algorithm"] == "")[0].size == self.NUM_ALGOS/2 else 1
        out.has_validation_o = 0 if np.where(stats_O["algorithm"] == "")[0].size == self.NUM_ALGOS else 1
        out.gage_type = gage_type.upper()
        if FLPEno== False:    
            #flpe dimenstion
            a_dim_flpe = out.createDimension("num_algos_flpe", None)
            c_dim_flpe = out.createDimension("nchar_flpe", None)
            c_dim_gage = out.createDimension("nchar_gage", None)
            t_dim_flpe = out.createDimension("time_flpe", len(stats_flpe["t"]))
            t_v_flpe = out.createVariable("time_flpe", "i4", ("time_flpe",))
            t_v_flpe.units = "days since Jan 1 Year 1"       
            t_v_flpe[:] = stats_flpe["t"]
            
            a_v_flpe = out.createVariable("algorithm_flpe", 'S1', ("num_algos_flpe", "nchar_flpe"),)        
            a_v_flpe[:] = stringtochar(stats_flpe["algorithm"][0].astype("S16"))
           
            gid_v_flpe = out.createVariable("gageID_flpe", "S1", ("num_algos_flpe", "nchar_gage"), fill_value=fill)
            gid_v_flpe[:] = stringtochar(stats_flpe["Gid"][:].astype("S16"))
            r_v_flpe = out.createVariable("Spearmanr_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            r_v_flpe[:] = np.where(np.isclose(stats_flpe["Spearmanr"], empty), fill, stats_flpe["Spearmanr"])
            sige_v_flpe = out.createVariable("SIGe_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            sige_v_flpe[:] = np.where(np.isclose(stats_flpe["SIGe"], empty), fill, stats_flpe["SIGe"])
            nse_v_flpe = out.createVariable("NSE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            nse_v_flpe[:] = np.where(np.isclose(stats_flpe["NSE"], empty), fill, stats_flpe["NSE"])
            rsq_v_flpe = out.createVariable("Rsq_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            rsq_v_flpe[:] = np.where(np.isclose(stats_flpe["Rsq"], empty), fill, stats_flpe["Rsq"])       
            kge_v_flpe = out.createVariable("KGE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            kge_v_flpe[:] = np.where(np.isclose(stats_flpe["KGE"], empty), fill, stats_flpe["KGE"])
            rmse_v_flpe = out.createVariable("RMSE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            rmse_v_flpe.units = "m^3/s"
            rmse_v_flpe[:] = np.where(np.isclose(stats_flpe["RMSE"], empty), fill, stats_flpe["RMSE"])
            n_v_flpe = out.createVariable("testn_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            n_v_flpe[:] = np.where(np.isclose(stats_flpe["n"], empty), fill, stats_flpe["n"])
            nrmse_v_flpe = out.createVariable("nRMSE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            nrmse_v_flpe.units = "none"
            nrmse_v_flpe[:] = np.where(np.isclose(stats_flpe["nRMSE"], empty), fill, stats_flpe["nRMSE"])
            nb_v_flpe = out.createVariable("nBIAS_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            nb_v_flpe.units = "none"
            nb_v_flpe[:] = np.where(np.isclose(stats_flpe["nBIAS"], empty), fill, stats_flpe["nBIAS"])
        else:
            #flpe fill dimenstion
            f_dim_flpe = out.createDimension("num_algos_flpe", 1)
            c_dim_flpe = out.createDimension("nchar_flpe", 1)
            c_dim_gage = out.createDimension("nchar_gage", 1)
            t_dim_flpe = out.createDimension("time_flpe", 1)
            t_v_flpe = out.createVariable("time_flpe", "i4", ("time_flpe",))
            t_v_flpe.units = "days since Jan 1 Year 1"       
            t_v_flpe[:] = empty
            
            a_v_flpe = out.createVariable("algorithm_flpe", 'S1', ("num_algos_flpe", "nchar_flpe"),)        
            a_v_flpe[:] = empty
           
            gid_v_flpe = out.createVariable("gageID_flpe", "S1", ("num_algos_flpe", "nchar_gage"), fill_value=fill)
            gid_v_flpe[:] = empty
            r_v_flpe = out.createVariable("Spearmanr_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            r_v_flpe[:] = empty
            sige_v_flpe = out.createVariable("SIGe_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            sige_v_flpe[:] = empty
            nse_v_flpe = out.createVariable("NSE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            nse_v_flpe[:] =empty
            rsq_v_flpe = out.createVariable("Rsq_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            rsq_v_flpe[:] = empty       
            kge_v_flpe = out.createVariable("KGE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            kge_v_flpe[:] = empty
            rmse_v_flpe = out.createVariable("RMSE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            rmse_v_flpe.units = "m^3/s"
            rmse_v_flpe[:] = empty
            n_v_flpe = out.createVariable("testn_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            n_v_flpe[:] = empty
            nrmse_v_flpe = out.createVariable("nRMSE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            nrmse_v_flpe.units = "none"
            nrmse_v_flpe[:] = empty
            nb_v_flpe = out.createVariable("nBIAS_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            nb_v_flpe.units = "none"
            nb_v_flpe[:] = empty
       
       
       
        if MOIno== False:
            #moi dimenstion
            a_dim_moi = out.createDimension("num_algos_moi", None)
            c_dim_moi = out.createDimension("nchar_moi", None)
            t_dim_moi = out.createDimension("time_moi", len(stats_moi["t"]))
            t_v_moi = out.createVariable("time_moi", "i4", ("time_moi",))
            t_v_moi.units = "days since Jan 1 Year 1"
            t_v_moi[:] = stats_flpe["t"]
            a_v_moi = out.createVariable("algorithm_moi", 'S1', ("num_algos_moi", "nchar_moi"),)
            a_v_moi[:] = stringtochar(stats_moi["algorithm"][0].astype("S16"))
            gid_v_moi = out.createVariable("gageID_moi", "S1", ("num_algos_moi", "nchar_gage"), fill_value=fill)
            gid_v_moi[:] = stringtochar(stats_moi["Gid"][:].astype("S16"))
            r_v_moi = out.createVariable("Spearmanr_moi", "f8", ("num_algos_moi",), fill_value=fill)
            r_v_moi[:] = np.where(np.isclose(stats_moi["Spearmanr"], empty), fill, stats_moi["Spearmanr"])            
            sige_v_moi = out.createVariable("SIGe_moi", "f8", ("num_algos_moi",), fill_value=fill)
            sige_v_moi[:] = np.where(np.isclose(stats_moi["SIGe"], empty), fill, stats_moi["SIGe"])        
            nse_v_moi = out.createVariable("NSE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            nse_v_moi[:] = np.where(np.isclose(stats_moi["NSE"], empty), fill, stats_moi["NSE"])       
            rsq_v_moi = out.createVariable("Rsq_moi", "f8", ("num_algos_moi",), fill_value=fill)
            rsq_v_moi[:] = np.where(np.isclose(stats_moi["Rsq"], empty), fill, stats_moi["Rsq"])                       
            kge_v_moi = out.createVariable("KGE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            kge_v_moi[:] = np.where(np.isclose(stats_moi["KGE"], empty), fill, stats_moi["KGE"])
            rmse_v_moi = out.createVariable("RMSE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            rmse_v_moi.units = "m^3/s"
            rmse_v_moi[:] = np.where(np.isclose(stats_moi["RMSE"], empty), fill, stats_moi["RMSE"])
            n_v_moi = out.createVariable("testn_moi", "f8", ("num_algos_moi",), fill_value=fill)
            n_v_moi[:] = np.where(np.isclose(stats_moi["n"], empty), fill, stats_moi["n"])
            nrmse_v_moi = out.createVariable("nRMSE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            nrmse_v_moi.units = "none"
            nrmse_v_moi[:] = np.where(np.isclose(stats_moi["nRMSE"], empty), fill, stats_moi["nRMSE"])
            nb_v_moi = out.createVariable("nBIAS_moi", "f8", ("num_algos_moi",), fill_value=fill)
            nb_v_moi.units = "none"
            nb_v_moi[:] = np.where(np.isclose(stats_moi["nBIAS"], empty), fill, stats_moi["nBIAS"])
        else:
            #moi dimenstion
            a_dim_moi = out.createDimension("num_algos_moi",1)
            c_dim_moi = out.createDimension("nchar_moi", 1)
            t_dim_moi = out.createDimension("time_moi", 1)
            t_v_moi = out.createVariable("time_moi", "i4", ("time_moi",))
            t_v_moi.units = "days since Jan 1 Year 1"
            t_v_moi[:] = empty
            a_v_moi = out.createVariable("algorithm_moi", 'S1', ("num_algos_moi", "nchar_moi"),)
            a_v_moi[:] =  empty
            gid_v_moi = out.createVariable("gageID_moi", "S1", ("num_algos_moi", "nchar_gage"), fill_value=fill)
            gid_v_moi[:] =  empty
            r_v_moi = out.createVariable("Spearmanr_moi", "f8", ("num_algos_moi",), fill_value=fill)
            r_v_moi[:] =  empty          
            sige_v_moi = out.createVariable("SIGe_moi", "f8", ("num_algos_moi",), fill_value=fill)
            sige_v_moi[:] =  empty       
            nse_v_moi = out.createVariable("NSE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            nse_v_moi[:] =  empty    
            rsq_v_moi = out.createVariable("Rsq_moi", "f8", ("num_algos_moi",), fill_value=fill)
            rsq_v_moi[:] = empty                       
            kge_v_moi = out.createVariable("KGE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            kge_v_moi[:] = empty
            rmse_v_moi = out.createVariable("RMSE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            rmse_v_moi.units = "m^3/s"
            rmse_v_moi[:] =  empty
            n_v_moi = out.createVariable("testn_moi", "f8", ("num_algos_moi",), fill_value=fill)
            n_v_moi[:] = np.where(np.isclose(stats_moi["n"], empty), fill, stats_moi["n"])
            nrmse_v_moi = out.createVariable("nRMSE_moi", "f8", ("num_algos_moi",), fill_value=fill)
            nrmse_v_moi.units = "none"
            nrmse_v_moi[:] =  empty
            nb_v_moi = out.createVariable("nBIAS_moi", "f8", ("num_algos_moi",), fill_value=fill)
            nb_v_moi.units = "none"
            nb_v_moi[:] =  empty 
       
     
        
         
        if OFFno== False:
        #offline dimenstion
            a_dim_o = out.createDimension("num_algos_o", None)
            c_dim_o = out.createDimension("nchar_o", None)
            t_dim_o = out.createDimension("time_o", len(stats_O["t"]))
            t_v_o = out.createVariable("time_o", "i4", ("time_o",))
            t_v_o.units = "days since Jan 1 Year 1"
            t_v_o[:] = stats_O["t"]
            a_v_o = out.createVariable("algorithm_o", 'S1', ("num_algos_o", "nchar_o"),)      
            a_v_o[:] = stringtochar(stats_O["algorithm"][0].astype("S16"))
            gid_v_o = out.createVariable("gageID_o", "S1", ("num_algos_o", "nchar_gage"), fill_value=fill)
            gid_v_o[:] = stringtochar(stats_O["Gid"][:].astype("S16"))
            r_v_o = out.createVariable("Spearmanr_o", "f8", ("num_algos_o",), fill_value=fill)
            r_v_o[:] = np.where(np.isclose(stats_O["Spearmanr"], empty), fill, stats_O["Spearmanr"])           
            sige_v_o = out.createVariable("SIGe_o", "f8", ("num_algos_o",), fill_value=fill)
            sige_v_o[:] = np.where(np.isclose(stats_O["SIGe"], empty), fill, stats_O["SIGe"])              
            nse_v_o = out.createVariable("NSE_o", "f8", ("num_algos_o",), fill_value=fill)
            nse_v_o[:] = np.where(np.isclose(stats_O["NSE"], empty), fill, stats_O["NSE"])           
            rsq_v_o = out.createVariable("Rsq_o", "f8", ("num_algos_o",), fill_value=fill)
            rsq_v_o[:] = np.where(np.isclose(stats_O["Rsq"], empty), fill, stats_O["Rsq"])       
            kge_v_o = out.createVariable("KGE_o", "f8", ("num_algos_o",), fill_value=fill)
            kge_v_o[:] = np.where(np.isclose(stats_O["KGE"], empty), fill, stats_O["KGE"])
            n_v_o = out.createVariable("testn_o", "f8", ("num_algos_o",), fill_value=fill)
            n_v_o[:] = np.where(np.isclose(stats_O["n"], empty), fill, stats_O["n"])
            rmse_v_o = out.createVariable("RMSE_o", "f8", ("num_algos_o",), fill_value=fill)
            rmse_v_o.units = "m^3/s"
            rmse_v_o[:] = np.where(np.isclose(stats_O["RMSE"], empty), fill, stats_O["RMSE"])
            n_v_o = out.createVariable("testn_o", "f8", ("num_algos_o",), fill_value=fill)
            n_v_o[:] = np.where(np.isclose(stats_O["n"], empty), fill, stats_O["n"])
            nrmse_v_o = out.createVariable("nRMSE_o", "f8", ("num_algos_o",), fill_value=fill)
            nrmse_v_o.units = "none"
            nrmse_v_o[:] = np.where(np.isclose(stats_O["nRMSE"], empty), fill, stats_O["nRMSE"])
            nb_v_o = out.createVariable("nBIAS_o", "f8", ("num_algos_o",), fill_value=fill)
            nb_v_o.units = "none"
            nb_v_o[:] = np.where(np.isclose(stats_O["nBIAS"], empty), fill, stats_O["nBIAS"])
        else:
             #offline dimenstion
            a_dim_o = out.createDimension("num_algos_o", 1)
            c_dim_o = out.createDimension("nchar_o", 1)
            t_dim_o = out.createDimension("time_o", 1)
            t_v_o = out.createVariable("time_o", "i4", ("time_o",))
            t_v_o.units = "days since Jan 1 Year 1"
            t_v_o[:] = empty
            a_v_o = out.createVariable("algorithm_o", 'S1', ("num_algos_o", "nchar_o"),)      
            a_v_o[:] = empty
            gid_v_o = out.createVariable("gageID_o", "S1", ("num_algos_o", "nchar_gage"), fill_value=fill)
            gid_v_o[:] = empty
            r_v_o = out.createVariable("Spearmanr_o", "f8", ("num_algos_o",), fill_value=fill)
            r_v_o[:] = empty           
            sige_v_o = out.createVariable("SIGe_o", "f8", ("num_algos_o",), fill_value=fill)
            sige_v_o[:] = empty               
            nse_v_o = out.createVariable("NSE_o", "f8", ("num_algos_o",), fill_value=fill)
            nse_v_o[:] = empty           
            rsq_v_o = out.createVariable("Rsq_o", "f8", ("num_algos_o",), fill_value=fill)
            rsq_v_o[:] = empty    
            kge_v_o = out.createVariable("KGE_o", "f8", ("num_algos_o",), fill_value=fill)
            kge_v_o[:] = empty           
            rmse_v_o = out.createVariable("RMSE_o", "f8", ("num_algos_o",), fill_value=fill)
            rmse_v_o.units = "m^3/s"
            rmse_v_o[:] =empty
            n_v_o = out.createVariable("testn_o", "f8", ("num_algos_o",), fill_value=fill)
            n_v_o[:] =empty
            nrmse_v_o = out.createVariable("nRMSE_o", "f8", ("num_algos_o",), fill_value=fill)
            nrmse_v_o.units = "none"
            nrmse_v_o[:] =empty
            nb_v_o = out.createVariable("nBIAS_o", "f8", ("num_algos_o",), fill_value=fill)
            nb_v_o.units = "none"
            nb_v_o[:] =empty        

       

        out.close()

   

def get_reach_data(input_json,index_to_run):
        """Retrun dictionary of reach data.
        
        Parameters
        ----------
        input_json: str
            string name of json file used to detect what to execute on
            
        Returns
        -------
        dictionary of reach data
        """
        
        if index_to_run == -235:
            index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
        else: 
            index=index_to_run

        with open(INPUT / input_json) as json_file:
            reach_data = json.load(json_file)[index]
        return reach_data

def run_validation():
    """Orchestrate validation operations."""

    try:
        reach_json = sys.argv[1]
        run_type = sys.argv[2]
    except IndexError:
        reach_json = "reaches.json"
        run_type = "unconstrained"

    try:
        index_to_run=int(sys.argv[3]) #integer
    except IndexError:
        index_to_run=-235 #code to indicate AWS run

    reach_data = get_reach_data(reach_json,index_to_run)

    vc = ValidationConfluence(reach_data,FLPE,MOI, OFFLINE, INPUT, OUTPUT, run_type)
    vc.validate()

if __name__ == "__main__": 

    start = datetime.datetime.now()
    run_validation()
    end = datetime.datetime.now()
    print(f"Execution time: {end - start}")
