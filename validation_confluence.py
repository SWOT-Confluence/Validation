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

# Local imports
from val.validation import stats

# Third-party imports
from netCDF4 import Dataset, stringtochar
import numpy as np

# Constants
INPUT = Path("/mnt/data/input")
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
    read_time_data()
        read time of observations from SWOT files
    validate()
        run validation operations on gage data and FLPE data; write stats
    write(stats, time, reach_id, gage_type)
        write stats to NetCDF file
    """

    INT_FILL = -999
    NUM_ALGOS = 10

    def __init__(self, reach_data, offline_dir, input_dir, output_dir, run_type):

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
        self.output_dir = output_dir


    def read_gage_data(self, sos_file):
        """Read gage data from SoS file and stores in gage data dictionary."""

        sos = Dataset(sos_file, 'r')

        # could be optimized using the global gage_agency varaible
        groups = list(sos.groups.keys())
        gage_data = {}
        if "usgs" in groups:
            gage_data = self.get_gage_q(sos, "usgs")
        elif "grdc" in groups:
            gage_data = self.get_gage_q(sos, "grdc")
        elif "Hidroweb" in groups:
            gage_data = self.get_gage_q(sos, "Hidroweb")
        elif "ABOM" in groups:
            gage_data = self.get_gage_q(sos, "ABOM")
        elif "WSC" in groups:
            gage_data = self.get_gage_q(sos, "WSC")
        elif "MLIT" in groups:
            gage_data = self.get_gage_q(sos, "MLIT")
        elif "DEFRA" in groups:
            gage_data = self.get_gage_q(sos, "DEFRA")
            if gage_data == {}:
                gage_data = self.get_gage_q(sos, "EAU")
            
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
        


            

        gage_data = {}
        if len(index[0]) != 0:
            if self.run_type == "constrained":
                # if constraind check and see if the gage selected at this index is a 0
                if gage["CAL"][:][index] == 1:
                    warnings.warn('gauge found was calibration.. This is an unconstrained run, so it will not be used for validation')
                    return gage_data
            gage_data["type"] = gage_type
            gage_data["q"] = gage[f"{gage_type}_q"][index][:].filled(np.nan)
            if gage_type == "usgs":
                gage_data["qt"] = gage[f"{gage_type}_qt"][index][:].filled(self.INT_FILL).astype(int)
            else:
                qt = gage[f"{gage_type}_qt"][index][:].astype(int) + 366    # Vt in matlab tatetime +366 converts to Python ordinal
                gage_data["qt"] = qt.filled(self.INT_FILL)
        return gage_data        

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
            "metro_q_uc":"dschg_m",
            "bam_q_uc":"dschg_b",
            "hivdi_q_uc":"dschg_h",
            "momma_q_uc":"dschg_o",
            "sads_q_uc":"dschg_s",
            # "sic4dvar_q_uc":"dschg_i",
            "consensus_q_uc":"dschg_c",
            "d_x_area":"d_x_area",
            "d_x_area_u":"d_x_area_u",
        }

        offline_file = f"{offline_dir}/{self.reach_id}_offline.nc"
        off = Dataset(offline_file, 'r')
        offline_data = {}
        offline_data["geobam_q_c"] = off[convention_dict["bam_q_c"]][:].filled(np.nan)
        offline_data["hivdi_q_c"] = off[convention_dict["hivdi_q_c"]][:].filled(np.nan)
        offline_data["metroman_q_c"] = off[convention_dict["metro_q_c"]][:].filled(np.nan)
        offline_data["momma_q_c"] = off[convention_dict["momma_q_c"]][:].filled(np.nan)
        offline_data["sad_q_c"] = off[convention_dict["sads_q_c"]][:].filled(np.nan)
        offline_data["geobam_q_uc"] = off[convention_dict["bam_q_uc"]][:].filled(np.nan)
        offline_data["hivdi_q_uc"] = off[convention_dict["hivdi_q_uc"]][:].filled(np.nan)
        offline_data["metroman_q_uc"] = off[convention_dict["metro_q_uc"]][:].filled(np.nan)
        offline_data["momma_q_uc"] = off[convention_dict["momma_q_uc"]][:].filled(np.nan)
        offline_data["sad_q_uc"] = off[convention_dict["sads_q_uc"]][:].filled(np.nan)
        offline_data["consensus_q_c"] = off[convention_dict["consensus_q_c"]][:].filled(np.nan)
        offline_data["consensus_q_uc"] = off[convention_dict["consensus_q_uc"]][:].filled(np.nan)
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
                print(time)
                warnings.warn('problem with time conversion to ordinal, most likely nan value')
        # return [ (epoch + datetime.timedelta(seconds=t)).toordinal() for t in time ]   # Check if this format works
        return ordinal_times

    def validate(self):
        """Run validation operations on gage data and FLPE data; write stats."""
        
        # SWOT time 
        time = self.read_time_data()

        # Data fill values
        data = {
            "algorithm": np.full((self.NUM_ALGOS), fill_value=""),
            "NSE": np.full((self.NUM_ALGOS), fill_value=-9999),
            "Rsq": np.full((self.NUM_ALGOS), fill_value=-9999),
            "KGE": np.full((self.NUM_ALGOS), fill_value=-9999),
            "RMSE": np.full((self.NUM_ALGOS), fill_value=-9999),
            "n": np.full((self.NUM_ALGOS), fill_value=-9999),
            "nRMSE":np.full((self.NUM_ALGOS), fill_value=-9999),
            "nBIAS":np.full((self.NUM_ALGOS), fill_value=-9999),
            "rRMSE":np.full((self.NUM_ALGOS), fill_value=-9999),
        }
        no_offline = False
        # Check if there is data to validate
        if self.gage_data:
            if self.offline_data:
                #### Check should go here for all nan gauge data ---------------------------------
                data = stats(time, self.offline_data, self.gage_data["qt"], 
                            self.gage_data["q"], str(self.reach_id), 
                            self.output_dir / "figs")
            else:
                warnings.warn('No offline data found...')
                no_offline = True
        else:
            warnings.warn('No gauge found for reach...')
            
        # Write out valid or invalid data
        gage_type = "No data" if not self.gage_data else self.gage_data["type"]

        if (gage_type != "No data") and (no_offline != True):
            self.write(data, self.reach_id, gage_type)

    def write(self, stats, reach_id, gage_type):
        """Write stats to NetCDF file.
        
        Parameters
        ----------
        stats: dict
            dictionary of stats for each algorithm
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
        out.has_validation = 0 if np.where(stats["algorithm"] == "")[0].size == self.NUM_ALGOS else 1
        out.gage_type = gage_type.upper()

        a_dim = out.createDimension("num_algos", None)
        c_dim = out.createDimension("nchar", None)
        t_dim = out.createDimension("time", len(stats["t"]))
        t_v = out.createVariable("time", "i4", ("time",))
        t_v.units = "days since Jan 1 Year 1"
        t_v[:] = stats["t"]

        a_v = out.createVariable("algorithm", 'S1', ("num_algos", "nchar"),)
        a_v[:] = stringtochar(stats["algorithm"].astype("S16"))
        
        nse_v = out.createVariable("NSE", "f8", ("num_algos",), fill_value=fill)
        nse_v[:] = np.where(np.isclose(stats["NSE"], empty), fill, stats["NSE"])

        rsq_v = out.createVariable("Rsq", "f8", ("num_algos",), fill_value=fill)
        rsq_v[:] = np.where(np.isclose(stats["Rsq"], empty), fill, stats["Rsq"])

        kge_v = out.createVariable("KGE", "f8", ("num_algos",), fill_value=fill)
        kge_v[:] = np.where(np.isclose(stats["KGE"], empty), fill, stats["KGE"])

        rmse_v = out.createVariable("RMSE", "f8", ("num_algos",), fill_value=fill)
        rmse_v.units = "m^3/s"
        rmse_v[:] = np.where(np.isclose(stats["RMSE"], empty), fill, stats["RMSE"])
        
        nrmse_v = out.createVariable("nRMSE", "f8", ("num_algos",), fill_value=fill)
        nrmse_v.units = "none"
        nrmse_v[:] = np.where(np.isclose(stats["nRMSE"], empty), fill, stats["nRMSE"])
        
        nb_v = out.createVariable("nBIAS", "f8", ("num_algos",), fill_value=fill)
        nb_v.units = "none"
        nb_v[:] = np.where(np.isclose(stats["nBIAS"], empty), fill, stats["nBIAS"])
        
        rrmse_v = out.createVariable("rRMSE", "f8", ("num_algos",), fill_value=fill)
        rrmse_v.units = "none"
        rrmse_v[:] = np.where(np.isclose(stats["rRMSE"], empty), fill, stats["rRMSE"])

        n_v = out.createVariable("testn", "f8", ("num_algos",), fill_value=fill)
        n_v[:] = np.where(np.isclose(stats["n"], empty), fill, stats["n"])

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

    vc = ValidationConfluence(reach_data, OFFLINE, INPUT, OUTPUT, run_type)
    vc.validate()

if __name__ == "__main__": 

    start = datetime.datetime.now()
    run_validation()
    end = datetime.datetime.now()
    print(f"Execution time: {end - start}")
