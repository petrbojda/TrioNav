#!/usr/bin/env python

import sys
sys.path.append("../../srcpy")

import numpy as np
from datacontainers import data_containers as dc
from utils import nav_plots as nplt
#import cProfile, pstats, io

# Path to the configuration file, where data specs are stored
data_config_file = "./data_specs.cnf"


def main(path_to_data):

    # Load Data from .mat file
    print ("Path to data:",path_to_data)

    IMU_measurements = dc.List_MP_IMU()
    IMU_measurements.append_from_m_file(data_path=path_to_data)
    print ("Number of measured points",IMU_measurements.__len__())
    print ("Counters interval",IMU_measurements.get_count_interval())
    print ("Duration of measurement",IMU_measurements.get_time_duration() /60,"min")

    sel = {"rotX_tp": None, "rotY_tp": None, "rotZ_tp": None, "DrotX_tp": None, "DrotY_tp": None, "DrotZ_tp": None ,
           "accX_tp": None, "accY_tp": None, "accZ_tp": None, "DaccX_tp": None, "DaccY_tp": None, "DaccZ_tp": None ,
           "time_tp":None, "cnt_tp": None}
    nplt.static_plot_IMU(IMU_measurements,sel,None)

    if IMU_measurements:

        Plot_data = IMU_measurements.get_array_data_sel(selection = sel)

    a = np.max(Plot_data["rotX"])
    print (a)




#pr = cProfile.Profile()
#pr.enable()

if __name__ == "__main__":
    path_to_data = dc.cnf_file_scenario_select(data_config_file)

if path_to_data:
        main(path_to_data)

#pr.disable()
#s = io.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print (s.getvalue())
