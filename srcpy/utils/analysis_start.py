#!/usr/bin/env python

import sys

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

    IMU_measurements = dcIMU.List_MP_IMU()
    IMU_measurements.append_from_m_file(data_path=path_to_data)
    print(60 * '-')
    print('IMU Measurements')
    print(60 * '-')
    print("Number of measured points", IMU_measurements.__len__())
    print("Counters interval", IMU_measurements.get_count_interval())
    print("SYStime interval", IMU_measurements.get_time_interval(),' in minutes: ', IMU_measurements.get_time_interval()[0]/60,IMU_measurements.get_time_interval()[1]/60 )
    print("Duration of measurement", IMU_measurements.get_time_duration() / 60, "min")
    print(60 * '-')

    selIMU = {  "rotX_tp": None, "rotY_tp": None, "rotZ_tp": None, "DrotX_tp": None, "DrotY_tp": None, "DrotZ_tp": None,
                "accX_tp": None, "accY_tp": None, "accZ_tp": None, "DaccX_tp": None, "DaccY_tp": None, "DaccZ_tp": None,
                "time_tp":None, "cnt_tp": None}

    pltimu.three_rotacc_raw(IMU_measurements,selIMU,None,', raw measurements',1)

    GPS_measurements = dcGPS.List_MP_GPS()
    GPS_measurements.append_from_m_file(data_path=path_to_data)
    print(60 * '-')
    print('GPS Measurements')
    print(60 * '-')
    print("Number of measured points", GPS_measurements.__len__())
    print("Counters interval", GPS_measurements.get_count_interval())
    print("SYStime interval", GPS_measurements.get_time_interval(),' in minutes: ', GPS_measurements.get_time_interval()[0]/60,GPS_measurements.get_time_interval()[1]/60 )
    print("Duration of measurement", GPS_measurements.get_time_duration() / 60, "min")
    print(60 * '-')

    selGPS = {  "lat_tp": None, "lon_tp": None, "alt_tp": None, "gr_spd_tp": None, "heading_tp": None,
                "dopp_vel_N_tp": None, "dopp_vel_E_tp": None, "dopp_vel_D_tp": None, "n_o_satellites_tp": None,
                "fix_val_tp": None,  "vdop_tp": None,  "hdop_tp": None,  "pdop_tp": None,
                "SYStime_tp":None, "GPStime_tp":None, "cnt_tp": None}

    pltgps.all_raw(GPS_measurements, selGPS, None,', raw measurements',3)
    
    RTKuBlox_measurements = dcRTK.List_MP_RTK()
    RTKuBlox_measurements.append_from_m_file(data_path=path_to_data, row_offset=31)
    print(60 * '-')
    print('RTKuBlox Measurements')
    print(60 * '-')
    print("Number of measured points", RTKuBlox_measurements.__len__())
    print("Counters interval", RTKuBlox_measurements.get_count_interval())
    print("SYStime interval", RTKuBlox_measurements.get_time_interval(),' in minutes: ', RTKuBlox_measurements.get_time_interval()[0]/60,RTKuBlox_measurements.get_time_interval()[1]/60 )
    print("Duration of measurement", RTKuBlox_measurements.get_time_duration() / 60, "min")
    print(60 * '-')

    selRTKuBlox = {  "lat_tp": None, "lon_tp": None, "alt_tp": None, "RTKage_tp": None, "stdUN_tp": None,
                "stdNorth_tp": None, "stdEast_tp": None, "stdUp_tp": None, "n_o_satellites_tp": None,
                "fix_val_tp": None,  "stdNE_tp": None,  "stdEU_tp": None,  "RTKratio_tp": None,
                "SYStime_tp":None, "GPStime_tp":None, "cnt_tp": None}

    pltrtk.all_raw(RTKuBlox_measurements, selRTKuBlox, None,', RTK uBlox raw',7)
    
    RTKGeoX91_measurements = dcRTK.List_MP_RTK()
    RTKGeoX91_measurements.append_from_m_file(data_path=path_to_data, row_offset=45)
    print(60 * '-')
    print('RTKGeoX91 Measurements')
    print(60 * '-')
    print("Number of measured points", RTKGeoX91_measurements.__len__())
    print("Counters interval", RTKGeoX91_measurements.get_count_interval())
    print("SYStime interval", RTKGeoX91_measurements.get_time_interval(),' in minutes: ', RTKGeoX91_measurements.get_time_interval()[0]/60,RTKGeoX91_measurements.get_time_interval()[1]/60 )
    print("Duration of measurement", RTKGeoX91_measurements.get_time_duration() / 60, "min")
    print(60 * '-')

    selRTKGeoX91 = {  "lat_tp": None, "lon_tp": None, "alt_tp": None, "RTKage_tp": None, "stdUN_tp": None,
                "stdNorth_tp": None, "stdEast_tp": None, "stdUp_tp": None, "n_o_satellites_tp": None,
                "fix_val_tp": None,  "stdNE_tp": None,  "stdEU_tp": None,  "RTKratio_tp": None,
                "SYStime_tp":None, "GPStime_tp":None, "cnt_tp": None}

    pltrtk.all_raw(RTKGeoX91_measurements, selRTKGeoX91, None,', RTK GeoX91 raw',11)




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
