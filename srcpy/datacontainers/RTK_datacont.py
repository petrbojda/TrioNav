import scipy.io as sio
import numpy as np



class MeasuredPoint_RTK(object):
    def __init__(self, SYStime=None, GPStime=None, cnt=None, lat=None, lon=None, alt=None, n_o_satellites=None, fix_val=None,
                 stdNorth=None, stdEast=None, stdUp=None, stdNE=None, stdEU=None, stdUN=None, RTKage=None, RTKratio=None):

        self._SYStime = SYStime
        self._GPStime = GPStime
        self._cnt = cnt
        self._lat = lat
        self._lon =lon
        self._alt = alt
        self._stdNorth =stdNorth
        self._stdEast = stdEast
        self._stdUp =stdUp
        self._stdNE = stdNE
        self._stdEU =stdEU
        self._stdUN = stdUN
        self._RTKage =RTKage
        self._RTKratio = RTKratio
        self._n_o_satellites =n_o_satellites
        self._fix_val =fix_val



class List_MP_RTK(list):

    def __init__(self,IMU_sampling_rate = 200):
        super().__init__()
        self._sampling_rate = IMU_sampling_rate
        
        self._cnt_interval = (0,0)
        self._SYStime_interval =  (0,0)
        self._GPStime_interval =  (0,0)
        
        self._lat_interval =  (0,0)
        self._lon_interval = (0,0)
        self._alt_interval = (0,0)
        
        self._stdNorth_interval = (0,0)
        self._stdEast_interval = (0,0)
        self._stdUp_interval = (0,0)
        self._stdNE_interval = (0,0)
        self._stdEU_interval = (0,0)
        self._stdUN_interval = (0,0)
        self._RTKage_interval = (0,0)
        self._RTKratio_interval = (0,0)

        self._n_o_satellites_interval = (0,0)
        self._fix_val_interval = (0,0)


    def append_from_m_file(self, data_path,row_offset = 31):
        measured_data_mfile = sio.loadmat(data_path)
        raw_data = measured_data_mfile["data"]
        no_d = len(raw_data)

        for itr in range(0, no_d ):
            #print('For itr=',itr,'the GPS time is',raw_data[itr, row_offset],'which converts to',raw_data[itr, row_offset] != 'NaN','since it is of type ', type(raw_data[itr, row_offset]))
            if ~np.isnan(raw_data[itr, row_offset]):
                self.append(MeasuredPoint_RTK ( cnt = int(raw_data[itr, 0]),
                                                SYStime = int(raw_data[itr, 0]) / self._sampling_rate,
                                                GPStime = int(raw_data[itr, row_offset+0]),
                                                lat=float(raw_data[itr, row_offset+1]),
                                                lon=float(raw_data[itr, row_offset+2]),
                                                alt=float(raw_data[itr, row_offset+3]),
                                                stdNorth=float(raw_data[itr, row_offset+6]),
                                                stdEast=float(raw_data[itr, row_offset+7]),
                                                stdUp=float(raw_data[itr, row_offset+8]),
                                                stdNE=float(raw_data[itr, row_offset+9]),
                                                stdEU=float(raw_data[itr, row_offset+10]),
                                                stdUN=float(raw_data[itr, row_offset+11]),
                                                RTKage=float(raw_data[itr, row_offset+12]),
                                                RTKratio=float(raw_data[itr, row_offset+13]),
                                                n_o_satellites=float(raw_data[itr, row_offset+5]),
                                                fix_val=float(raw_data[itr, row_offset+4])))

        if raw_data[itr, row_offset+0]:
            self._cnt_interval = (min([elem._cnt for elem in self]),max([elem._cnt for elem in self]))
            self._SYStime_interval = (min([elem._SYStime for elem in self]),max([elem._SYStime for elem in self]))
            self._GPStime_interval = (min([elem._GPStime for elem in self]),max([elem._GPStime for elem in self]))

            self._lat_interval = (min([elem._lat for elem in self]),max([elem._lat for elem in self]))
            self._lon_interval = (min([elem._lon for elem in self]),max([elem._lon for elem in self]))
            self._alt_interval = (min([elem._alt for elem in self]),max([elem._alt for elem in self]))
            self._stdNorth_interval = (min([elem._stdNorth for elem in self]),max([elem._stdNorth for elem in self]))
            self._stdEast_interval = (min([elem._stdEast for elem in self]),max([elem._stdEast for elem in self]))
            self._stdUp_interval = (min([elem._stdUp for elem in self]),max([elem._stdUp for elem in self]))
            self._stdNE_interval = (min([elem._stdNE for elem in self]),max([elem._stdNE for elem in self]))

            self._stdEU_interval = (min([elem._stdEU for elem in self]),max([elem._stdEU for elem in self]))
            self._stdUN_interval = (min([elem._RTKage for elem in self]),max([elem._RTKage for elem in self]))
            self._RTKage_interval = (min([elem._stdUN for elem in self]),max([elem._stdUN for elem in self]))
            self._RTKratio_interval = (min([elem._RTKratio for elem in self]),max([elem._RTKratio for elem in self]))

            self._n_o_satellites_interval = (min([elem._n_o_satellites for elem in self]),max([elem._n_o_satellites for elem in self]))
            self._fix_val_interval = (min([elem._fix_val for elem in self]),max([elem._fix_val for elem in self]))


    def get_count_interval(self):
        self._cnt_interval = (min([elem._cnt for elem in self]),max([elem._cnt for elem in self]))
        return self._cnt_interval

    def get_time_interval(self):
        self._SYStime_interval = (min([elem._SYStime for elem in self]),max([elem._SYStime for elem in self]))
        return self._SYStime_interval

    def get_time_duration(self):
        time_duration = (self._SYStime_interval[1] - self._SYStime_interval[0])
        return time_duration

    def get_array_data_sel(self, **kwarg):

        if 'cnt' in kwarg:
            cnt_i = kwarg['cnt'] if (len(kwarg['cnt']) == 2) else (kwarg['cnt'],kwarg['cnt'])
        else:
            cnt_i = self._cnt_interval

        if 'SYStime' in kwarg:
            SYStime_i = kwarg['SYStime'] if (len(kwarg['SYStime']) == 2) else (kwarg['SYStime'],kwarg['SYStime'])
        else:
            SYStime_i = self._SYStime_interval

        if 'GPStime' in kwarg:
            GPStime_i = kwarg['GPStime'] if (len(kwarg['GPStime']) == 2) else (kwarg['GPStime'],kwarg['GPStime'])
        else:
            GPStime_i = self._GPStime_interval

        if 'lat' in kwarg:
            lat_i = kwarg['lat'] if (len(kwarg['lat']) == 2) else (kwarg['lat'],kwarg['lat'])
        else:
            lat_i = self._lat_interval

        if 'lon' in kwarg:
            lon_i = kwarg['lon'] if (len(kwarg['lon']) == 2) else (kwarg['lon'],kwarg['lon'])
        else:
            lon_i = self._lon_interval

        if 'alt' in kwarg:
            alt_i = kwarg['alt'] if (len(kwarg['alt']) == 2) else (kwarg['alt'],kwarg['alt'])
        else:
            alt_i = self._alt_interval

        if 'stdNorth' in kwarg:
            stdNorth_i = kwarg['stdNorth'] if (len(kwarg['stdNorth']) == 2) else (kwarg['stdNorth'],kwarg['stdNorth'])
        else:
            stdNorth_i = self._stdNorth_interval

        if 'stdEast' in kwarg:
            stdEast_i = kwarg['stdEast'] if (len(kwarg['stdEast']) == 2) else (kwarg['stdEast'],kwarg['stdEast'])
        else:
            stdEast_i = self._stdEast_interval

        if 'stdUp' in kwarg:
            stdUp_i = kwarg['stdUp'] if (len(kwarg['stdUp']) == 2) else (kwarg['stdUp'],kwarg['stdUp'])
        else:
            stdUp_i = self._stdUp_interval

        if 'stdNE' in kwarg:
            stdNE_i = kwarg['stdNE'] if (len(kwarg['stdNE']) == 2) else (kwarg['stdNE'],kwarg['stdNE'])
        else:
            stdNE_i = self._stdNE_interval

        if 'stdEU' in kwarg:
            stdEU_i = kwarg['stdEU'] if (len(kwarg['stdEU']) == 2) else (kwarg['stdEU'],kwarg['stdEU'])
        else:
            stdEU_i = self._stdEU_interval

        if 'stdUN' in kwarg:
            stdUN_i = kwarg['stdUN'] if (len(kwarg['stdUN']) == 2) else (kwarg['stdUN'],kwarg['stdUN'])
        else:
            stdUN_i = self._stdUN_interval

        if 'RTKage' in kwarg:
            RTKage_i = kwarg['RTKage'] if (len(kwarg['RTKage']) == 2) else (kwarg['RTKage'],kwarg['RTKage'])
        else:
            RTKage_i = self._RTKage_interval

        if 'RTKratio' in kwarg:
            RTKratio_i = kwarg['RTKratio'] if (len(kwarg['RTKratio']) == 2) else (kwarg['RTKratio'],kwarg['RTKratio'])
        else:
            RTKratio_i = self._RTKratio_interval

        if 'n_o_satellites' in kwarg:
            n_o_satellites_i = kwarg['n_o_satellites'] if (len(kwarg['n_o_satellites']) == 2) else (kwarg['n_o_satellites'],kwarg['n_o_satellites'])
        else:
            n_o_satellites_i = self._n_o_satellites_interval

        if 'fix_val' in kwarg:
            fix_val_i = kwarg['fix_val'] if (len(kwarg['fix_val']) == 2) else (kwarg['fix_val'],kwarg['fix_val'])
        else:
            fix_val_i = self._fix_val_interval




        if 'selection' in kwarg:
            cnt_i = kwarg['selection']['cnt_tp'] if kwarg['selection']['cnt_tp'] else self._cnt_interval
            SYStime_i = kwarg['selection']['SYStime_tp'] if kwarg['selection']['SYStime_tp'] else self._SYStime_interval
            GPStime_i = kwarg['selection']['GPStime_tp'] if kwarg['selection']['GPStime_tp'] else self._GPStime_interval

            fix_val_i = kwarg['selection']['fix_val_tp'] if kwarg['selection']['fix_val_tp'] else self._fix_val_interval
            n_o_satellites_i = kwarg['selection']['n_o_satellites_tp'] if kwarg['selection']['n_o_satellites_tp'] else self._n_o_satellites_interval

            lat_i = kwarg['selection']['lat_tp'] if kwarg['selection']['lat_tp'] else self._lat_interval
            lon_i = kwarg['selection']['lon_tp'] if kwarg['selection']['lon_tp'] else self._lon_interval
            alt_i = kwarg['selection']['alt_tp'] if kwarg['selection']['alt_tp'] else self._alt_interval
            stdNorth_i = kwarg['selection']['stdNorth_tp'] if kwarg['selection']['stdNorth_tp'] else self._stdNorth_interval

            stdEU_i = kwarg['selection']['stdEU_tp'] if kwarg['selection']['stdEU_tp'] else self._stdEU_interval
            stdUN_i = kwarg['selection']['stdUN_tp'] if kwarg['selection']['stdUN_tp'] else self._stdUN_interval
            RTKage_i = kwarg['selection']['RTKage_tp'] if kwarg['selection']['RTKage_tp'] else self._RTKage_interval
            RTKratio_i = kwarg['selection']['RTKratio_tp'] if kwarg['selection']['RTKratio_tp'] else self._RTKratio_interval

            stdNE_i = kwarg['selection']['stdNE_tp'] if kwarg['selection']['stdNE_tp'] else self._stdNE_interval
            stdEast_i = kwarg['selection']['stdEast_tp'] if kwarg['selection']['stdEast_tp'] else self._stdEast_interval
            stdUp_i = kwarg['selection']['stdUp_tp'] if kwarg['selection']['stdUp_tp'] else self._stdUp_interval

        cnt_sel = [elem._cnt for elem in self if         (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        SYStime_sel = [elem._SYStime for elem in self if (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        GPStime_sel = [elem._GPStime for elem in self if (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        fix_val_sel = [elem._fix_val for elem in self if (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        n_o_satellites_sel = [elem._n_o_satellites for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        lat_sel = [elem._lat for elem in self if        (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        lon_sel = [elem._lon for elem in self if        (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        alt_sel = [elem._alt for elem in self if        (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        stdNorth_sel = [elem._stdNorth for elem in self if( SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        stdEU_sel = [elem._stdEU for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        stdUN_sel = [elem._stdUN for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        RTKage_sel = [elem._RTKage for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        RTKratio_sel = [elem._RTKratio for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        stdUp_sel = [elem._stdUp for elem in self if    (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        stdEast_sel = [elem._stdEast for elem in self if (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        stdNE_sel = [elem._stdNE for elem in self if      ( SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            stdNorth_i[0] <= elem._stdNorth <= stdNorth_i[1] and
                                                            stdEU_i[0] <= elem._stdEU <= stdEU_i[1] and
                                                            stdUN_i[0] <= elem._stdUN <= stdUN_i[1] and
                                                            RTKage_i[0] <= elem._RTKage <= RTKage_i[1] and
                                                            RTKratio_i[0] <= elem._RTKratio <= RTKratio_i[1] and
                                                            stdUp_i[0] <= elem._stdUp <= stdUp_i[1] and
                                                            stdNE_i[0] <= elem._stdNE <= stdNE_i[1] and
                                                            stdEast_i[0] <= elem._stdEast <= stdEast_i[1])]

        RTK_data = {  "SYStime": np.array(SYStime_sel),
                      "GPStime": np.array(GPStime_sel),
                      "cnt": np.array(cnt_sel),
                      "fix_val": np.array(fix_val_sel),
                      "n_o_satellites": np.array(n_o_satellites_sel),
                      "lat": np.array(lat_sel),
                      "lon": np.array(lon_sel),
                      "alt": np.array(alt_sel),
                      "stdNorth": np.array(stdNorth_sel),
                      "stdEU": np.array(stdEU_sel),
                      "stdUN": np.array(stdUN_sel),
                      "RTKage": np.array(RTKage_sel),
                      "RTKratio": np.array(RTKratio_sel),
                      "stdEast": np.array(stdEast_sel),
                      "stdUp": np.array(stdUp_sel),
                      "stdNE": np.array(stdNE_sel)}

        return RTK_data

