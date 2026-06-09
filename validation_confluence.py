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
import argparse
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
from sos_read.sos_read import download_sos

# Third-party imports
from netCDF4 import Dataset, stringtochar, chartostring
import numpy as np

# Constants
INPUT = Path("/mnt/data/input")
FLPE = Path("/mnt/data/flpe")
MOI = Path("/mnt/data/moi")
OFFLINE = Path("/mnt/data/offline")
OUTPUT = Path("/mnt/data/output")
TMP_DIR = Path("/tmp")

# SVS validation paths (set to None to use SoS priors instead)
SVS_FILE = Path("/mnt/data/input/svs/SVS_v1_0.5_trans_v17.nc")
SVS_EXCLUDE_JSON = Path("/mnt/data/input/svs/reachids_svs_gages_used_to_train.json")

FLPE_MOI_ALGOS = [
    "metroman",
    "busboi",
    "hivdi",
    "momma",
    "sad",
    "sic4dvar",
    "consensus",
]

MOI_BASE_ALGOS = [
    "metroman",
    "busboi",
    "hivdi",
    "momma",
    "sad",
    "sic4dvar",
]


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
    NUM_ALGOS = len(FLPE_MOI_ALGOS)  # flpe/moi: metroman, busboi, hivdi, momma, sad, sic4dvar, consensus
    NUM_ALGOS_OFFLINE = 16

    def __init__(self, reach_data, flpe_dir, moi_dir, offline_dir, input_dir, output_dir, run_type, gage_dir, svs_file=None, exclude_json=None):
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
        gage_dir: Path
            path to priors SOS directory
        svs_file: str or Path, optional
            path to SVS NetCDF file. If provided, SVS gage data is used
            instead of SoS gage data.
        exclude_json: str or Path, optional
            path to JSON file with reach IDs to exclude from SVS validation
        """
        
        self.input_dir = input_dir
        self.run_type = run_type
        self.reach_id = reach_data["reach_id"]
        print('Processing', self.reach_id)
        if svs_file is not None:
            self.gage_data = self.read_gage_data_svs(svs_file, exclude_json)
        else:
            self.gage_data = self.read_gage_data(gage_dir / reach_data["sos"])

        #turn off offline for this run (v4)
        self.offline_data = {}
        
        self.flpe_data = self.read_flpe_data(flpe_dir)
        try:
            self.moi_data = self.read_moi_data(moi_dir)
        except FileNotFoundError:
            warnings.warn(f'No MOI file found for reach {self.reach_id}, skipping MOI validation')
            self.moi_data = {}
        self.output_dir = output_dir

    def read_gage_data(self, sos_file):
        """Read gage data from SoS file and stores in gage data dictionary."""

        sos = Dataset(sos_file, 'r')
        gage_data = {}
        for gage_agency in sos.gauge_agency.split(';'):
            gage_data = self.get_gage_q(sos, gage_agency)
            if gage_data != {}:
                print('found gage')
                break

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
        print('here is index we are working with ', index)
        # if its more than one, we take it down to a scalar
        if len(index[0]) > 1:
            warnings.warn('multiple gages for this reach. Selecting closest meanQ to model')
            #pull model q for this reach
            modelindex = np.where(self.reach_id == sos['reaches']['reach_id'][:].filled(np.nan))
            model_q = sos['model']['mean_q'][modelindex][:].filled(np.nan)
            gmq = []
            glt = []
            for Gindex in index:
                #pull mean q and timeseries lenghts
                gmq.append(gage[f"{gage_type}_mean_q"][Gindex][:].filled(np.nan))
                t = gage[f"{gage_type}_qt"][Gindex][:].filled(self.INT_FILL).astype(int)
                glt.append(len(t[t > 0]))
                 
            if np.isnan(model_q):
                #when model is nan, choose longest timeseries
                index = np.array(index[np.argmax(np.array(glt))])
                if np.size(index) > 1:
                    warnings.warn('model was nan and times are same length')
                    index = index[0]
            else:
                #othewise closest mean
                index = np.array(index[np.argmin(np.abs(np.array(glt) - model_q))])
                if np.size(index) > 1:
                    warnings.warn('identical mean q values')
                    index = index[0]
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

    def read_gage_data_svs(self, svs_file, exclude_json=None):
        """Read gage data from SVS NetCDF file and return gage data dictionary.
        
        Same output format as read_gage_data: a dictionary with keys
        'type', 'q', 'qt', 'gid', or an empty dict if no match found.
        
        Parameters
        ----------
        svs_file: str or Path
            path to the SVS NetCDF file
        exclude_json: str or Path, optional
            path to JSON file containing reach IDs used for ML training
            that should be excluded from validation
            
        Returns
        -------
        dict
            gage data dictionary matching read_gage_data output format
        """
        
        svs = Dataset(svs_file, 'r')
        gage_data = self.get_gage_q_svs(svs, exclude_json)
        svs.close()
        
        if gage_data:
            print('found SVS gage')
        else:
            print('no SVS gage match for reach', self.reach_id)
        
        return gage_data

    def get_gage_q_svs(self, svs, exclude_json=None):
        """Return discharge and discharge time from SVS NetCDF dataset.
        
        Matches self.reach_id against reach_id_v17, excludes training IDs,
        and returns the same dictionary format as get_gage_q.
        
        Parameters
        ----------
        svs: NetCDF dataset
            SVS NetCDF dataset reference
        exclude_json: str or Path, optional
            path to JSON file containing reach IDs to exclude
        
        Returns
        -------
        dict
            dictionary with keys 'type', 'q', 'qt', 'gid' or empty dict
        """
        
        # Load exclude list if provided
        exclude_ids = set()
        if exclude_json is not None:
            try:
                with open(exclude_json, 'r') as f:
                    data = json.load(f)
                if isinstance(data, dict):
                    exclude_ids = set(int(v) for v in data.values())
                elif isinstance(data, list):
                    exclude_ids = set(int(v) for v in data)
            except Exception as e:
                warnings.warn(f'Could not load exclude JSON: {e}')
        
        # Check if this reach should be excluded
        if int(self.reach_id) in exclude_ids:
            warnings.warn(f'Reach {self.reach_id} is in the training exclusion list, skipping')
            return {}
        
        # Read reach IDs - shape is (station, num_rchs)
        svs_reaches = svs['reach_id_v17'][:]
        
        if svs_reaches.ndim == 2:
            # Search across all num_rchs columns for a match
            svs_rid = svs_reaches.filled(np.nan)
            # Find station indices where any column matches self.reach_id
            match_mask = np.any(svs_rid == self.reach_id, axis=1)
            station_indices = np.where(match_mask)[0]
        else:
            svs_rid = svs_reaches.filled(np.nan)
            station_indices = np.where(svs_rid == self.reach_id)[0]
        
        print('SVS station indices for reach:', station_indices)
        
        if len(station_indices) == 0:
            return {}
        
        # If multiple stations match, pick the one with the longest valid Q record
        if len(station_indices) > 1:
            warnings.warn(f'Multiple SVS stations match reach {self.reach_id}, selecting longest record')
            best_idx = None
            best_count = -1
            for si in station_indices:
                q_series = svs['Q'][si, :].filled(np.nan)
                valid_count = np.count_nonzero(~np.isnan(q_series) & (q_series > 0))
                if valid_count > best_count:
                    best_count = valid_count
                    best_idx = si
            station_idx = best_idx
        else:
            station_idx = station_indices[0]
        
        # Extract Q for this station - shape (time,)
        q = svs['Q'][station_idx, :].filled(np.nan)
        
        # Convert SVS time to ordinal days
        # SVS time variable is in days; date_ymd gives [year, month, day] for each timestep
        date_ymd = svs['date_ymd'][:, :]  # shape (3, time) — [ymd_component, time]
        
        qt = np.full(len(q), self.INT_FILL, dtype=int)
        for i in range(date_ymd.shape[1]):
            try:
                year = int(date_ymd[0, i])
                month = int(date_ymd[1, i])
                day = int(date_ymd[2, i])
                qt[i] = datetime.date(year, month, day).toordinal()
            except:
                qt[i] = self.INT_FILL
        
        # Extract station ID as string
        # station_id is shape (strlen_station, station, string1)
        try:
            sid_chars = svs['station_id'][:, station_idx, 0]
            if hasattr(sid_chars, 'filled'):
                sid_chars = sid_chars.filled(b' ')
            gid = b''.join(sid_chars).decode('utf-8', errors='ignore').strip()
        except:
            gid = str(svs['station'][station_idx])
        
        gage_data = {}
        gage_data["type"] = "SVS"
        gage_data["q"] = q
        gage_data["qt"] = qt
        gage_data["gid"] = gid
        
        return gage_data
    
    def read_moi_data(self, moi_dir):
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

        def safe_read_q(group_name):
            if group_name in moi.groups and "q" in moi[group_name].variables:
                return moi[f"{group_name}/q"][:].filled(np.nan)
            return -9999

        for algo in MOI_BASE_ALGOS:
            moi_data[algo] = safe_read_q(algo)

        moi.close()

        # MOI output does not write a consensus group, so compute it here
        # from the available basin-scale algorithm discharge series.
        ref_key = next((k for k in MOI_BASE_ALGOS if not np.isscalar(moi_data[k])), None)
        if ref_key is not None:
            allq = np.full((len(MOI_BASE_ALGOS), len(moi_data[ref_key])), np.nan)
            for row, algo in enumerate(MOI_BASE_ALGOS):
                algv = moi_data[algo]
                if np.isscalar(algv) and algv == -9999:
                    continue
                algv = algv.copy()
                algv[algv < 0] = np.nan
                allq[row, :] = algv
            moi_data["consensus"] = np.nanmedian(allq, axis=0)
        else:
            moi_data["consensus"] = -9999

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
            if np.isscalar(v) and v == -9999:
                invalid += 1
            elif not np.isscalar(v) and np.count_nonzero(~np.isnan(v)) == 0:
                invalid += 1
        if invalid == self.NUM_ALGOS:
            print('moi IS NOT VALID')
            return False
        return True
        
    def read_flpe_data(self, flpe_dir):
        """Reads data from flpe module and returns dictionary.
        
        Parameters
        ----------
        flpe_dir: Path
            path to flpe data directory
        
        Returns
        -------
        dictionary of algorithm flpe results
        """
        convention_dict = {
            "metroman": "average/allq",
            "busboi": "q/q",
            "hivdi": "reach/Q",
            "momma": "Q",
            "sad": "Qa",
            "sic4dvar": "Q_da",
            "consensus": "consensus_q",
        }

        flpe_file_map = {
            "metroman": f"{flpe_dir}/metroman/{self.reach_id}_metroman.nc",
            "busboi": f"{flpe_dir}/busboi/{self.reach_id}_busboi.nc",
            "hivdi": f"{flpe_dir}/hivdi/{self.reach_id}_hivdi.nc",
            "momma": f"{flpe_dir}/momma/{self.reach_id}_momma.nc",
            "sad": f"{flpe_dir}/sad/{self.reach_id}_sad.nc",
            "sic4dvar": f"{flpe_dir}/sic4dvar/{self.reach_id}_sic4dvar.nc",
            "consensus": f"{flpe_dir}/consensus/{self.reach_id}_consensus.nc",
        }

        flpe_data = {}
        conlen = 0

        for algo in FLPE_MOI_ALGOS:
            flpe_file = flpe_file_map[algo]
            try:
                flpe_ds = Dataset(flpe_file, 'r')
            except Exception:
                flpe_data[algo] = -9999
                continue

            flpe_data[algo] = flpe_ds[convention_dict[algo]][:].filled(np.nan)
            conlen = len(flpe_data[algo])
            flpe_ds.close()

        if conlen > 0:
            if self.is_flpe_valid(flpe_data):
                return flpe_data
            else:
                return {}
        else:
            return {}

    def is_flpe_valid(self, flpe_data):
        """Check if flpe data is only comprised of NaN values.
        
        Returns
        -------
        False if all NaN values are present otherwise True
        """
        
        invalid = 0
        for v in flpe_data.values():
            if np.isscalar(v) and v == -9999:
                invalid += 1
            elif not np.isscalar(v) and np.count_nonzero(~np.isnan(v)) == 0:
                invalid += 1
        if invalid == self.NUM_ALGOS:
            print('flpe IS NOT VALID')
            return False
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
            "metro_q_c": "dschg_gm",
            "bam_q_c": "dschg_gb",
            "boi_q_c": "dschg_ga",
            "hivdi_q_c": "dschg_gh",
            "momma_q_c": "dschg_go",
            "sads_q_c": "dschg_gs",
            "sic4dvar_q_c": "dschg_gi",
            "consensus_q_c": "dschg_gc",
            "metro_q_uc": "dschg_m",
            "sic4dvar_q_uc": "dschg_i",
            "bam_q_uc": "dschg_b",
            "boi_q_uc": "dschg_a",
            "hivdi_q_uc": "dschg_h",
            "momma_q_uc": "dschg_o",
            "sads_q_uc": "dschg_s",
            "consensus_q_uc": "dschg_c",
            "d_x_area": "d_x_area",
            "d_x_area_u": "d_x_area_u",
        }

        offline_file = f"{offline_dir}/{self.reach_id}_offline.nc"
        off = Dataset(offline_file, 'r')
        offline_data = {}
        offline_data[convention_dict["bam_q_c"]] = off[convention_dict["bam_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["boi_q_c"]] = off[convention_dict["boi_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["hivdi_q_c"]] = off[convention_dict["hivdi_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["metro_q_c"]] = off[convention_dict["metro_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["momma_q_c"]] = off[convention_dict["momma_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["sads_q_c"]] = off[convention_dict["sads_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["sic4dvar_q_c"]] = off[convention_dict["sic4dvar_q_c"]][:].filled(np.nan)
        offline_data[convention_dict["sic4dvar_q_uc"]] = off[convention_dict["sic4dvar_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["bam_q_uc"]] = off[convention_dict["bam_q_uc"]][:].filled(np.nan)
        offline_data[convention_dict["boi_q_uc"]] = off[convention_dict["boi_q_uc"]][:].filled(np.nan)
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
            if np.count_nonzero(~np.isnan(v)) == 0:
                invalid += 1
        if invalid == self.NUM_ALGOS_OFFLINE:
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
        epoch = datetime.datetime(2000, 1, 1, 0, 0, 0)
        ordinal_times = []

        for t in time:
            try:
                ordinal_times.append((epoch + datetime.timedelta(seconds=t)).toordinal())
            except:
                ordinal_times.append(np.nan)
                warnings.warn('problem with time conversion to ordinal, most likely nan value')
        return ordinal_times

    def validate(self):
        """Run validation operations on gage data and FLPE data; write stats."""
        # SWOT time 
        time = self.read_time_data()
        algo_dim = int(self.NUM_ALGOS)
        Tdim = len(time)
        # Data fill values
        no_offline = True
        data_O = {
            "algorithm": np.full((self.NUM_ALGOS_OFFLINE), fill_value=""),
            "Gid": np.full((self.NUM_ALGOS_OFFLINE), fill_value=""),
            "pearsonr": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "SIGe": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "NSE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "Rsq": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "KGE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "RMSE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "n": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "nRMSE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "nBIAS": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "t": np.full(Tdim, fill_value=-9999),
            "consensus": np.full(Tdim, fill_value=-9999),
        }

        no_flpe = False
        # Check if there is data to validate
        if self.gage_data:
            try:
                if self.flpe_data:
                    data_flpe = stats(time, self.flpe_data, self.gage_data["qt"], 
                                      self.gage_data["q"], self.gage_data["gid"], str(self.reach_id), 
                                      self.output_dir / "figs")
                else:
                    warnings.warn('No flpe data found...')
                    no_flpe = True
            except Exception as e:
                warnings.warn(f'stats() failed for flpe reach {self.reach_id}: {e}')
                no_flpe = True
        else:
            warnings.warn('No gauge found for reach...')

        data_moi = {
            "algorithm": np.full(algo_dim, fill_value=""),
            "Gid": np.full(algo_dim, fill_value=""),
            "pearsonr": np.full(algo_dim, fill_value=-9999),
            "SIGe": np.full(algo_dim, fill_value=-9999),
            "NSE": np.full(algo_dim, fill_value=-9999),           
            "Rsq": np.full(algo_dim, fill_value=-9999),
            "KGE": np.full(algo_dim, fill_value=-9999),           
            "RMSE": np.full(algo_dim, fill_value=-9999),           
            "n": np.full(algo_dim, fill_value=-9999),           
            "nRMSE": np.full(algo_dim, fill_value=-9999),           
            "nBIAS": np.full(algo_dim, fill_value=-9999),
            "t": np.full(Tdim, fill_value=-9999),
            "consensus": np.full(Tdim, fill_value=-9999),
        }

        no_moi = False
        # Check if there is data to validate
        if self.gage_data:
            if self.moi_data:
                data_moi = stats(time, self.moi_data, self.gage_data["qt"], 
                                 self.gage_data["q"], self.gage_data["gid"], str(self.reach_id), 
                                 self.output_dir / "figs")
            else:
                warnings.warn('No moi data found...')
                no_moi = True
        else:
            warnings.warn('No gauge found for reach...')

        data_O = {
            "algorithm": np.full((self.NUM_ALGOS_OFFLINE), fill_value=""),
            "Gid": np.full((self.NUM_ALGOS_OFFLINE), fill_value=""),
            "pearsonr": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "SIGe": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "NSE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),           
            "Rsq": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "KGE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),           
            "RMSE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),           
            "n": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),           
            "nRMSE": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),           
            "nBIAS": np.full((self.NUM_ALGOS_OFFLINE), fill_value=-9999),
            "t": np.full(Tdim, fill_value=-9999),
            "consensus": np.full(Tdim, fill_value=-9999),
        }
        
        no_offline = False
        # Check if there is data to validate
        if self.gage_data:
            if self.offline_data:
                data_O = stats(time, self.offline_data, self.gage_data["qt"], 
                               self.gage_data["q"], self.gage_data["gid"], str(self.reach_id), 
                               self.output_dir / "figs")
            else:
                warnings.warn('No offline data found...')
                no_offline = True
        else:
            warnings.warn('No gauge found for reach...')
            
        # Write out valid or invalid data
        gage_type = "No data" if not self.gage_data else self.gage_data["type"]       
        ALLnone = np.all([no_flpe, no_moi, no_offline])     
        if (gage_type != "No data") and (ALLnone != True):
            self.write(data_flpe, data_moi, data_O, self.reach_id, gage_type, [no_flpe, no_moi, no_offline])

    def write(self, stats_flpe, stats_moi, stats_O, reach_id, gage_type, GO):
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

        FLPEno = GO[0]
        MOIno = GO[1]
        OFFno = GO[2]

        fill = -999999999999
        empty = -9999
        
        def safe_gid(gid_array):
            """Flatten Gid array to 1D array of strings for stringtochar."""
            flat = []
            for g in gid_array:
                if isinstance(g, np.ndarray):
                    flat.append(str(g.flat[0]))
                elif isinstance(g, (list, tuple)):
                    flat.append(str(g[0]))
                else:
                    flat.append(str(g))
            return np.array(flat)
        
        out = Dataset(self.output_dir / "stats" / f"{reach_id}_validation.nc", 'w')
        out.reach_id = reach_id
        out.description = f"Statistics for reach: {reach_id}"
        out.history = datetime.datetime.now().strftime('%d-%b-%Y %H:%M:%S')
        out.has_validation_flpe = 0 if np.where(stats_flpe["algorithm"] == "")[0].size == self.NUM_ALGOS else 1
        out.has_validation_moi  = 0 if np.where(stats_moi["algorithm"]  == "")[0].size == self.NUM_ALGOS else 1
        out.has_validation_o    = 0 if np.where(stats_O["algorithm"]    == "")[0].size == self.NUM_ALGOS_OFFLINE else 1
        out.gage_type = gage_type.upper()
        
        # Separate fixed dimensions for flpe/moi and offline
        out.createDimension("num_algos_flpe", self.NUM_ALGOS)
        out.createDimension("num_algos_offline", self.NUM_ALGOS_OFFLINE)
        c_dim_flpe = out.createDimension("nchar_flpe", None)
        c_dim_gage = out.createDimension("nchar_gage", None)

        time_values = None
        for stats_dict in (stats_flpe, stats_moi, stats_O):
            candidate = np.asarray(stats_dict["t"])
            if candidate.size > 0 and not np.all(np.isclose(candidate, empty)):
                time_values = candidate
                break
        if time_values is None:
            time_values = np.array([], dtype=int)

        t_dim = out.createDimension("time", len(time_values))
        t_v_flpe = out.createVariable("time", "i4", ("time",))
        t_v_flpe.units = "days since Jan 1 Year 1"
        t_v_flpe[:] = time_values

        # --- FLPE variables (use num_algos_flpe) ---
        if FLPEno == False:    
            a_v_flpe = out.createVariable("algorithm_flpe", 'S1', ("num_algos_flpe", "nchar_flpe"),)        
            a_v_flpe[:] = stringtochar(stats_flpe["algorithm"][0].astype("S16"))           
            gid_v_flpe = out.createVariable("gageID_flpe", "S1", ("num_algos_flpe", "nchar_gage"), fill_value=fill)
            gid_v_flpe[:] = stringtochar(stats_flpe["Gid"][:].astype("S16"))
            r_v_flpe = out.createVariable("pearsonr_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            r_v_flpe[:] = np.where(np.isclose(stats_flpe["pearsonr"], empty), fill, stats_flpe["pearsonr"])
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
            consensus_flpe = out.createVariable("consensus_flpe", "f8", ("time",), fill_value=fill)
            consensus_flpe.units = "m^3/s"
            consensus_flpe[:] = np.where(np.isclose(stats_flpe["consensus"], empty), fill, stats_flpe["consensus"])
        else:
            a_v_flpe = out.createVariable("algorithm_flpe", 'S1', ("num_algos_flpe", "nchar_flpe"),)        
            a_v_flpe[:] = empty           
            gid_v_flpe = out.createVariable("gageID_flpe", "S1", ("num_algos_flpe", "nchar_gage"), fill_value=fill)
            gid_v_flpe[:] = empty
            r_v_flpe = out.createVariable("pearsonr_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            r_v_flpe[:] = empty
            sige_v_flpe = out.createVariable("SIGe_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            sige_v_flpe[:] = empty
            nse_v_flpe = out.createVariable("NSE_flpe", "f8", ("num_algos_flpe",), fill_value=fill)
            nse_v_flpe[:] = empty
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
            consensus_flpe = out.createVariable("consensus_flpe", "f8", ("time",), fill_value=fill)
            consensus_flpe.units = "m^3/s"
            consensus_flpe[:] = empty

        # --- MOI variables (use num_algos_flpe) ---
        if MOIno == False:
            a_v_moi = out.createVariable("algorithm_moi", 'S1', ("num_algos_flpe", "nchar_flpe"),)
            a_v_moi[:] = stringtochar(stats_moi["algorithm"][0].astype("S16"))
            gid_v_moi = out.createVariable("gageID_moi", "S1", ("num_algos_flpe", "nchar_gage"), fill_value=fill)
            gid_v_moi[:] = stringtochar(stats_moi["Gid"][:].astype("S16"))
            r_v_moi = out.createVariable("pearsonr_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            r_v_moi[:] = np.where(np.isclose(stats_moi["pearsonr"], empty), fill, stats_moi["pearsonr"])            
            sige_v_moi = out.createVariable("SIGe_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            sige_v_moi[:] = np.where(np.isclose(stats_moi["SIGe"], empty), fill, stats_moi["SIGe"])        
            nse_v_moi = out.createVariable("NSE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            nse_v_moi[:] = np.where(np.isclose(stats_moi["NSE"], empty), fill, stats_moi["NSE"])       
            rsq_v_moi = out.createVariable("Rsq_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            rsq_v_moi[:] = np.where(np.isclose(stats_moi["Rsq"], empty), fill, stats_moi["Rsq"])                       
            kge_v_moi = out.createVariable("KGE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            kge_v_moi[:] = np.where(np.isclose(stats_moi["KGE"], empty), fill, stats_moi["KGE"])
            rmse_v_moi = out.createVariable("RMSE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            rmse_v_moi.units = "m^3/s"
            rmse_v_moi[:] = np.where(np.isclose(stats_moi["RMSE"], empty), fill, stats_moi["RMSE"])
            n_v_moi = out.createVariable("testn_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            n_v_moi[:] = np.where(np.isclose(stats_moi["n"], empty), fill, stats_moi["n"])
            nrmse_v_moi = out.createVariable("nRMSE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            nrmse_v_moi.units = "none"
            nrmse_v_moi[:] = np.where(np.isclose(stats_moi["nRMSE"], empty), fill, stats_moi["nRMSE"])
            nb_v_moi = out.createVariable("nBIAS_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            nb_v_moi.units = "none"
            nb_v_moi[:] = np.where(np.isclose(stats_moi["nBIAS"], empty), fill, stats_moi["nBIAS"])
        else:
            a_v_moi = out.createVariable("algorithm_moi", 'S1', ("num_algos_flpe", "nchar_flpe"),)
            a_v_moi[:] = empty
            gid_v_moi = out.createVariable("gageID_moi", "S1", ("num_algos_flpe", "nchar_gage"), fill_value=fill)
            gid_v_moi[:] = empty
            r_v_moi = out.createVariable("pearsonr_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            r_v_moi[:] = empty          
            sige_v_moi = out.createVariable("SIGe_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            sige_v_moi[:] = empty       
            nse_v_moi = out.createVariable("NSE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            nse_v_moi[:] = empty    
            rsq_v_moi = out.createVariable("Rsq_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            rsq_v_moi[:] = empty                       
            kge_v_moi = out.createVariable("KGE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            kge_v_moi[:] = empty
            rmse_v_moi = out.createVariable("RMSE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            rmse_v_moi.units = "m^3/s"
            rmse_v_moi[:] = empty
            n_v_moi = out.createVariable("testn_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            n_v_moi[:] = empty
            nrmse_v_moi = out.createVariable("nRMSE_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            nrmse_v_moi.units = "none"
            nrmse_v_moi[:] = empty
            nb_v_moi = out.createVariable("nBIAS_moi", "f8", ("num_algos_flpe",), fill_value=fill)
            nb_v_moi.units = "none"
            nb_v_moi[:] = empty

        # --- Offline variables (use num_algos_offline) ---
        if OFFno == False:
            a_v_o = out.createVariable("algorithm_o", 'S1', ("num_algos_offline", "nchar_flpe"),)      
            a_v_o[:] = stringtochar(stats_O["algorithm"][0].astype("S16"))
            gid_v_o = out.createVariable("gageID_o", "S1", ("num_algos_offline", "nchar_gage"), fill_value=fill)
            gid_v_o[:] = stringtochar(stats_O["Gid"][:].astype("S16"))
            r_v_o = out.createVariable("pearsonr_o", "f8", ("num_algos_offline",), fill_value=fill)
            r_v_o[:] = np.where(np.isclose(stats_O["pearsonr"], empty), fill, stats_O["pearsonr"])                  
            sige_v_o = out.createVariable("SIGe_o", "f8", ("num_algos_offline",), fill_value=fill)
            sige_v_o[:] = np.where(np.isclose(stats_O["SIGe"], empty), fill, stats_O["SIGe"])                
            nse_v_o = out.createVariable("NSE_o", "f8", ("num_algos_offline",), fill_value=fill)
            nse_v_o[:] = np.where(np.isclose(stats_O["NSE"], empty), fill, stats_O["NSE"])           
            rsq_v_o = out.createVariable("Rsq_o", "f8", ("num_algos_offline",), fill_value=fill)
            rsq_v_o[:] = np.where(np.isclose(stats_O["Rsq"], empty), fill, stats_O["Rsq"])       
            kge_v_o = out.createVariable("KGE_o", "f8", ("num_algos_offline",), fill_value=fill)
            kge_v_o[:] = np.where(np.isclose(stats_O["KGE"], empty), fill, stats_O["KGE"])          
            rmse_v_o = out.createVariable("RMSE_o", "f8", ("num_algos_offline",), fill_value=fill)
            rmse_v_o.units = "m^3/s"
            rmse_v_o[:] = np.where(np.isclose(stats_O["RMSE"], empty), fill, stats_O["RMSE"])
            n_v_o = out.createVariable("testn_o", "f8", ("num_algos_offline",), fill_value=fill)
            n_v_o[:] = np.where(np.isclose(stats_O["n"], empty), fill, stats_O["n"])
            nrmse_v_o = out.createVariable("nRMSE_o", "f8", ("num_algos_offline",), fill_value=fill)
            nrmse_v_o.units = "none"
            nrmse_v_o[:] = np.where(np.isclose(stats_O["nRMSE"], empty), fill, stats_O["nRMSE"])
            nb_v_o = out.createVariable("nBIAS_o", "f8", ("num_algos_offline",), fill_value=fill)
            nb_v_o.units = "none"
            nb_v_o[:] = np.where(np.isclose(stats_O["nBIAS"], empty), fill, stats_O["nBIAS"])
        else:
            a_v_o = out.createVariable("algorithm_o", 'S1', ("num_algos_offline", "nchar_flpe"),)      
            a_v_o[:] = empty
            gid_v_o = out.createVariable("gageID_o", "S1", ("num_algos_offline", "nchar_gage"), fill_value=fill)
            gid_v_o[:] = empty
            r_v_o = out.createVariable("pearsonr_o", "f8", ("num_algos_offline",), fill_value=fill)
            r_v_o[:] = empty           
            sige_v_o = out.createVariable("SIGe_o", "f8", ("num_algos_offline",), fill_value=fill)
            sige_v_o[:] = empty               
            nse_v_o = out.createVariable("NSE_o", "f8", ("num_algos_offline",), fill_value=fill)
            nse_v_o[:] = empty           
            rsq_v_o = out.createVariable("Rsq_o", "f8", ("num_algos_offline",), fill_value=fill)
            rsq_v_o[:] = empty    
            kge_v_o = out.createVariable("KGE_o", "f8", ("num_algos_offline",), fill_value=fill)
            kge_v_o[:] = empty           
            rmse_v_o = out.createVariable("RMSE_o", "f8", ("num_algos_offline",), fill_value=fill)
            rmse_v_o.units = "m^3/s"
            rmse_v_o[:] = empty
            n_v_o = out.createVariable("testn_o", "f8", ("num_algos_offline",), fill_value=fill)
            n_v_o[:] = empty
            nrmse_v_o = out.createVariable("nRMSE_o", "f8", ("num_algos_offline",), fill_value=fill)
            nrmse_v_o.units = "none"
            nrmse_v_o[:] = empty
            nb_v_o = out.createVariable("nBIAS_o", "f8", ("num_algos_offline",), fill_value=fill)
            nb_v_o.units = "none"
            nb_v_o[:] = empty   

        out.close()


def get_reach_data(input_json, index_to_run, sos_bucket):
    """Return dictionary of reach data.
    
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
        index = index_to_run

    with open(INPUT / input_json) as json_file:
        reach_data = json.load(json_file)[index]

    if sos_bucket:
        sos_file = TMP_DIR.joinpath(reach_data["sos"])
        download_sos(sos_bucket, sos_file)

    return reach_data

  
def create_args():
    """Create and return argparsers with command line arguments."""
    
    arg_parser = argparse.ArgumentParser(description='Integrate FLPE')
    arg_parser.add_argument('-i',
                            '--index',
                            type=int,
                            help='Index to specify input data to execute on')
    arg_parser.add_argument('-r',
                            '--reachjson',
                            type=str,
                            help='Name of the reaches.json',
                            default='reaches.json')
    arg_parser.add_argument('-t',
                            '--runtype',
                            type=str,
                            help='Indicates constrained or unconstrained run',
                            choices=['constrained', 'unconstrained'],
                            default='unconstrained')
    arg_parser.add_argument('-s',
                            '--sosbucket',
                            type=str,
                            help='Name of the SoS bucket and key to download from',
                            default='')

    return arg_parser


def run_validation():
    """Orchestrate validation operations."""
    
    # commandline arguments
    arg_parser = create_args()
    args = arg_parser.parse_args()

    reach_json = args.reachjson
    run_type = args.runtype
    sos_bucket = args.sosbucket

    index_to_run = args.index

    # 0.2 specify index to run. pull from command line arg or set to default = AWSf
    if args.index == -235:
        index_to_run = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))

    print('index_to_run: ', index_to_run)
    print('reach_json: ', reach_json)
    print('run_type: ', run_type)
    print('sos_bucket: ', sos_bucket)

    reach_data = get_reach_data(reach_json, index_to_run, sos_bucket)

    if sos_bucket:
        gage_dir = TMP_DIR
    else:
        gage_dir = INPUT.joinpath("sos")

    # Auto-detect SVS vs SoS: use SVS if the file exists, otherwise fall back to SoS priors
    svs_file = None
    exclude_json = None
    if SVS_FILE.exists():
        print(f'SVS file found: {SVS_FILE}. Remove SVS file in mnt/input/svs or change validation_confluence.py SVS filepaths to None to revert to sos priors.')
        svs_file = SVS_FILE
        if SVS_EXCLUDE_JSON.exists():
            exclude_json = SVS_EXCLUDE_JSON

    vc = ValidationConfluence(reach_data, FLPE, MOI, OFFLINE, INPUT, OUTPUT, run_type, gage_dir,
                              svs_file=svs_file, exclude_json=exclude_json)
    vc.validate()


if __name__ == "__main__": 

    start = datetime.datetime.now()
    run_validation()
    end = datetime.datetime.now()
    print(f"Execution time: {end - start}")
