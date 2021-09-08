"""Module to run validation operations and output stats.

Runs on a continent of data and requires JSON data for continent retrieved by
AWS Batch index.

TODO
- Temporarily only runs on Sacramento reaches -> need to remove this restriction

Class
-----
ValidationConfluence: Stores data and executes validation operations.

Constants
---------
INPUT_DIR: Path
    path to input directory
OUTPUT_DIR: Path
    path to output directory

Functions
---------
get_cont_data(cont_json, index)
    extract and return the continent data needs to be extracted for
run_validation()
    orchestrate validation operations
"""

# Standard imports
import json
import os
from os import scandir
from pathlib import Path
import sys

# Local imports
from validation import stats

# Third-party imports
from netCDF4 import Dataset, stringtochar
import numpy as np

# Constants
INPUT = Path("/home/nikki/Documents/confluence/workspace/validation/data/input/input")
OFFLINE_DIR = Path("/home/nikki/Documents/confluence/workspace/validation/data/input/offline")
OUTPUT = Path("/home/nikki/Documents/confluence/workspace/validation/data/output")
FIG_DIR = Path("/home/nikki/Documents/confluence/workspace/validation/data/figs")

class ValidationConfluence:
    """Class that runs validation operations for Confluence workflow.
    
    Attributes
    ----------
    cont: dict
        continent name key with associated numeric identifier values (list)
    fig_dir: Path
        path to figure directory (storage of hydrographs)
    gage_data: dict
        dictionary of gage reach identifiers, q, and qt
    input_dir: Path
        path to input directory
    NUM_ALGOS: int
        number of algorithms to store data for
    offline_ids: numpy.ndarray
        list of reach identifiers with offline data
    offline_data: dict
        dictionary of offline discharge values and time
    offline_dir: Path
        path to offline data directory
    output_dir: Path
        path to output directory
    reach_ids: numpy.ndarray
        array of reach identifiers for continent
    sos_file: Path
        path to SoS file for continent
    
    Methods
    -------
    get_reach_ids(sos_file)
        retreive reach identifiers from SoS file
    read_gage_data()
        read gage data from SoS file
    read_offline_data(reach_id)
        reads data from offline module and stores in flpe_data dictionary
    validate()
        run validation operations on gage data and FLPE data; write stats
    """

    NUM_ALGOS = 10

    def __init__(self, cont_json, offline_dir, sos_file, input_dir, output_dir,
        fig_dir):

        """
        Parameters
        ----------
        cont_json: dict
            continent name key with associated numeric identifier values (list)
        offline_dir: Path
            path to offline data directory
        sos_file: Path
            path to SoS file for continent
        input_dir: Path
            path to input directory
        output_dir: Path
            path to output directory
        fig_dir: Path
            path to figure directory (storage of hydrographs)
        """
        self.cont = cont_json
        self.fig_dir = fig_dir
        self.gage_data = {}
        self.input_dir = input_dir
        self.offline_ids = self.__get_offline_ids(offline_dir)
        self.offline_data = {}
        self.offline_dir = offline_dir
        self.output_dir = output_dir
        self.reach_ids = self.get_reach_ids(sos_file)
        self.sos_file = sos_file

    def __get_offline_ids(self, offline_dir):
        """Return reach identifiers that have offline data."""

        with scandir(offline_dir) as entries:
            reach_ids = [ int(entry.name.split('_')[0]) for entry in entries ]
        return(np.array(reach_ids))

    def get_reach_ids(self, sos_file):
        """Retreive reach identifiers from SoS file.
        
        Parameters
        ----------
        sos_file: Path
            path to continent SoS file

        Returns
        -------
        numpy.ndarray of reach identifieres
        """

        nc = Dataset(sos_file, 'r')
        rids = nc["reaches"]["reach_id"][:]
        nc.close()
        return rids

    def read_gage_data(self):
        """Read gage data from SoS file and stores in gage data dictionary."""

        sos = Dataset(self.sos_file, 'r')
        groups = list(sos["model"].groups.keys())
        
        if "usgs" in groups:
            usgs = sos["model"]["usgs"]
            self.gage_data["type"] = "usgs"
            self.gage_data["reach_id"] = usgs["usgs_reach_id"][:].filled(np.nan)
            self.gage_data["q"] = usgs["usgs_q"][:].filled(np.nan)
            self.gage_data["qt"] = usgs["usgs_qt"][:].filled(np.nan).astype(int)
        elif "grdc" in groups:
            grdc = sos["model"]["grdc"]
            self.gage_data["type"] = "grdc"
            self.gage_data["reach_id"] = grdc["grdc_reach_id"][:].filled(np.nan)
            self.gage_data["q"] = grdc["grdc_q"][:].filled(np.nan)
            self.gage_data["qt"] = grdc["grdc_qt"][:].filled(np.nan).astype(int) + 366    # Vt in matlab tatetime +366 converts to Python ordinal

        sos.close()

    def read_offline_data(self, reach_id):
        """Reads data from offline module and stores in flpe_data dictionary.
        
        Parameters
        ----------
        reach_id: int
            reach identifier for data file
        """

        offline_file = f"{self.offline_dir}/{reach_id}_offline.nc"
        off = Dataset(offline_file, 'r')
        self.offline_data["geobam_q_c"] = off["bam_q_c"][:]
        self.offline_data["hivdi_q_c"] = off["hivdi_q_c"][:]
        self.offline_data["metroman_q_c"] = off["metro_q_c"][:]
        self.offline_data["momma_q_c"] = off["momma_q_c"][:]
        self.offline_data["sad_q_c"] = off["sads_q_c"][:]
        self.offline_data["geobam_q_uc"] = off["bam_q_uc"][:]
        self.offline_data["hivdi_q_uc"] = off["hivdi_q_uc"][:]
        self.offline_data["metroman_q_uc"] = off["metro_q_uc"][:]
        self.offline_data["momma_q_uc"] = off["momma_q_uc"][:]
        self.offline_data["sad_q_uc"] = off["sads_q_uc"][:]
        off.close()

    def read_time_data(self, reach_id):
        """Read time of observations from SWOT files.

        Parameters
        ----------
        reach_id: int
            unique reach identifier
        
        Returns
        -------
        list of ordinal times
        """

        swot = Dataset(self.input_dir / "swot" / f"{reach_id}_SWOT.nc", 'r')
        time = swot["nt"][:].filled(np.nan)
        swot.close()
        return [ datetime.strptime(str(t), "%Y%m%d").toordinal() + 1 for t in time ]


    def validate(self, sac_rids):
        """Run validation operations on gage data and FLPE data; write stats.
        
        TODO
        - Remove restriction that only runs on sac reach identifiers
        """

        # for reach in self.reach_ids:    ## TODO
        for reach in sac_rids:
            data = {
                "algorithm": np.full((self.NUM_ALGOS), fill_value=""),
                "NSE": np.full((self.NUM_ALGOS), fill_value=-9999),
                "Rsq": np.full((self.NUM_ALGOS), fill_value=-9999),
                "KGE": np.full((self.NUM_ALGOS), fill_value=-9999),
                "RMSE": np.full((self.NUM_ALGOS), fill_value=-9999),
                "n": np.full((self.NUM_ALGOS), fill_value=-9999)
            }
            time = self.read_time_data(reach)
            # Check if gage data is present
            if reach in self.gage_data["reach_id"]:
                # Determine if offline data for gage data
                if reach in self.offline_ids:
                    # Gage data
                    index = np.where(self.gage_data["reach_id"] == reach)
                    qt = self.gage_data["qt"][index,:].flatten()
                    q = self.gage_data["q"][index,:].flatten()

                    # Offline data
                    self.read_offline_data(reach)
                    
                    # Stats
                    data = stats(time, self.offline_data, qt, q, str(reach),
                        self.fig_dir)
            self.write(data, time, reach, self.gage_data["type"])

    def write(self, stats, time, reach_id, gage_type):
        """Write stats to NetCDF file.
        
        Parameters
        ----------
        stats: dict
            dictionary of stats for each algorithm
        time: list
            list of observation times
        reach_id: int
            reach identifier for stats
        gage_type: str
            type of gage data used for validation
        """

        fill = -999999999999
        empty = -9999

        out = Dataset(self.output_dir / f"{reach_id}_validation.nc", 'w')
        out.reach_id = reach_id
        out.description = f"Statistics for reach: {reach_id}"
        out.history = datetime.now().strftime('%d-%b-%Y %H:%M:%S')
        out.has_validation = 0 if np.where(stats["algorithm"] == "")[0].size == self.NUM_ALGOS else 1
        out.gage_type = gage_type.upper()

        a_dim = out.createDimension("num_algos", None)
        c_dim = out.createDimension("nchar", None)
        t_dim = out.createDimension("time", len(time))
        t_v = out.createVariable("time", "f4", ("time",))
        t_v.units = "days since Jan 1 Year 1"
        t_v[:] = time

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

        n_v = out.createVariable("testn", "f8", ("num_algos",), fill_value=fill)
        n_v[:] = np.where(np.isclose(stats["n"], empty), fill, stats["n"])

        out.close()

def get_cont_data(cont_json, index):
    """Extract and return the continent data needs to be extracted for.
    
    Parameters
    ----------
    cont_json : str
        path to the file that contains the list of continents
    index: int
        integer index value to select json data
    
    Returns
    -------
    dictionary element of continent abbreviation and numeric identifiers
    """

    with open(cont_json) as jsonfile:
        data = json.load(jsonfile)
    return data[index]

def run_validation():
    """Orchestrate validation operations.
    
    TODO
    - Remove restriction of only running on sac data
    """

    # index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
    index = 3    ## TODO

    try:
        cont_json = sys.argv[1]
    except IndexError:
        cont_json = "continent.json"
    cont_data = get_cont_data(INPUT / cont_json, index)

    sos_file = INPUT / "sos" / f"{list(cont_data.keys())[0]}_apriori_rivers_v07_SOS.nc"
    vc = ValidationConfluence(cont_data, OFFLINE_DIR, sos_file, INPUT, OUTPUT, FIG_DIR)
    vc.read_gage_data()
    
    ## TODO sac reach ids
    with open(INPUT / "reaches.json") as json_file:
        data = json.load(json_file)
    sac_rids = [ e["reach_id"] for e in data ]
    
    vc.validate(sac_rids)

if __name__ == "__main__":

    from datetime import datetime    

    start = datetime.now()
    run_validation()
    end = datetime.now()
    print(f"Execution time: {end - start}")