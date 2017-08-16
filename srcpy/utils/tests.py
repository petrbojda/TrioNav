#!/usr/bin/env python

import sys
import types

sys.path.append("../../srcpy")

from datacontainers import IMU_datacont as dcIMU
from datacontainers import GPS_datacont as dcGPS
from datacontainers import RTK_datacont as dcRTK
from plots import IMU_plots as pltimu
from plots import GPS_plots as pltgps
from plots import RTK_plots as pltrtk
from utils import config_tools as cfgt

# import cProfile, pstats, io

# Path to the configuration file, where data specs are stored
data_config_file = "./data_specs.cnf"

def main(path_to_data):
    # Load Data from .mat file
    print("Path to data:", path_to_data)

    print('IMU Measurements')
    print(60 * '-')

    IMU_measurements = dcIMU.List_MP_IMU()
    print("Number of measured points before m-file reading: ", IMU_measurements.__len__())
    IMU_measurements.append_from_m_file(data_path=path_to_data)
    print(60 * '-')

    selIMU = {  "rotX_tp": None, "rotY_tp": None, "rotZ_tp": None, "DrotX_tp": None, "DrotY_tp": None, "DrotZ_tp": None,
                "accX_tp": None, "accY_tp": None, "accZ_tp": None, "DaccX_tp": None, "DaccY_tp": None, "DaccZ_tp": None,
                "time_tp":None, "cnt_tp": (0,1000)}

    IMU_selection = dcIMU.List_MP_IMU()
    IMUdata = IMU_measurements.get_array_data_sel(selection=selIMU)
    IMU_selection.append_from_dict(IMUdata)

    print("Number of measured points", IMU_measurements.__len__())
    print("Counters interval", IMU_measurements.get_count_interval())
    print("SYStime interval", IMU_measurements.get_time_interval(),' in minutes: ', IMU_measurements.get_time_interval()[0]/60,IMU_measurements.get_time_interval()[1]/60 )
    print("Duration of measurement", IMU_measurements.get_time_duration() / 60, "min")
    print(60 * '-')
    list_ = [ (n,v.__doc__) for n,v in IMU_measurements[1].__dict__.items() ]
    print('Items in a class: ',list_)
    print(60 * '-')
    print(60 * '-')
    print('IMU Selection')
    print(60 * '-')
    print(60 * '-')
    print("Number of selected points", IMU_selection.__len__())
    print("Counters interval", IMU_selection.get_count_interval())
    print("SYStime interval", IMU_selection.get_time_interval(),' in minutes: ', IMU_selection.get_time_interval()[0]/60,IMU_selection.get_time_interval()[1]/60 )
    print("Duration of measurement", IMU_selection.get_time_duration() / 60, "min")
    print(60 * '-')

    list_ = [ (n,v.__doc__) for n,v in IMU_selection[1].__dict__.items() ]
    print('Items in a class: ',list_)

    selIMU = {  "rotX_tp": None, "rotY_tp": None, "rotZ_tp": None, "DrotX_tp": None, "DrotY_tp": None, "DrotZ_tp": None,
                "accX_tp": None, "accY_tp": None, "accZ_tp": None, "DaccX_tp": None, "DaccY_tp": None, "DaccZ_tp": None,
                "time_tp":None, "cnt_tp": None}

    pltimu.three_rotacc_raw(IMU_selection,selIMU,None,', raw measurements',1)

# pr = cProfile.Profile()
# pr.enable()

if __name__ == "__main__":
    path_to_data = cfgt.cnf_file_scenario_select(data_config_file)

if path_to_data:
    main(path_to_data)

# pr.disable()
# s = io.StringIO()
# sortby = 'cumulative'
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
# print (s.getvalue())
