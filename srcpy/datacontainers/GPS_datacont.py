import scipy.io as sio
import numpy as np



class MeasuredPoint_GPS(object):
    def __init__(self, SYStime=None, GPStime=None, cnt=None, lat=None, lon=None, alt=None, gr_spd=None, vdop=None, hdop=None,
                 pdop=None, dopp_vel_D=None, dopp_vel_N=None, dopp_vel_E=None, heading=None, n_o_satellites=None, fix_val=None):

        self._SYStime = SYStime
        self._GPStime = GPStime
        self._cnt = cnt
        self._lat = lat
        self._lon =lon
        self._alt = alt
        self._gr_spd =gr_spd
        self._vdop = vdop
        self._hdop =hdop
        self._pdop = pdop
        self._dopp_vel_D =dopp_vel_D
        self._dopp_vel_N = dopp_vel_N
        self._dopp_vel_E =dopp_vel_E
        self._heading = heading
        self._n_o_satellites =n_o_satellites
        self._fix_val =fix_val



class List_MP_GPS(list):

    def __init__(self,IMU_sampling_rate = 200):
        super().__init__()
        self._sampling_rate = IMU_sampling_rate
        self._cnt_interval = (0,0)

        self._SYStime_interval =  (0,0)
        self._GPStime_interval =  (0,0)

        self._lat_interval =  (0,0)
        self._lon_interval = (0,0)
        self._alt_interval = (0,0)
        self._gr_spd_interval = (0,0)

        self._vdop_interval = (0,0)
        self._hdop_interval = (0,0)
        self._pdop_interval = (0,0)
        self._dopp_vel_N_interval = (0,0)
        self._dopp_vel_E_interval = (0,0)
        self._dopp_vel_D_interval = (0,0)
        self._heading_interval = (0,0)

        self._n_o_satellites_interval = (0,0)
        self._fix_val_interval = (0,0)


    def append_from_m_file(self, data_path,row_offset = 16):
        measured_data_mfile = sio.loadmat(data_path)
        raw_data = measured_data_mfile["data"]
        no_d = len(raw_data)

        for itr in range(0, no_d ):
            #print('For itr=',itr,'the GPS time is',raw_data[itr, row_offset],'which converts to',raw_data[itr, row_offset] != 'NaN','since it is of type ', type(raw_data[itr, row_offset]))
            if ~np.isnan(raw_data[itr, row_offset]):
                self.append(MeasuredPoint_GPS ( cnt = int(raw_data[itr, 0]),
                                                SYStime = int(raw_data[itr, 0]) / self._sampling_rate,
                                                GPStime = int(raw_data[itr, row_offset+0]),
                                                lat=float(raw_data[itr, row_offset+3]),
                                                lon=float(raw_data[itr, row_offset+4]),
                                                alt=float(raw_data[itr, row_offset+5]),
                                                gr_spd=float(raw_data[itr, row_offset+6]),
                                                vdop=float(raw_data[itr, row_offset+7]),
                                                hdop=float(raw_data[itr, row_offset+8]),
                                                pdop=float(raw_data[itr, row_offset+9]),
                                                dopp_vel_N=float(raw_data[itr, row_offset+11]),
                                                dopp_vel_E=float(raw_data[itr, row_offset+12]),
                                                dopp_vel_D=float(raw_data[itr, row_offset+10]),
                                                heading=float(raw_data[itr, row_offset+13]),
                                                n_o_satellites=float(raw_data[itr, row_offset+2]),
                                                fix_val=float(raw_data[itr, row_offset+1])))

        if raw_data[itr, row_offset+0]:
            self._cnt_interval = (min([elem._cnt for elem in self]),max([elem._cnt for elem in self]))
            self._SYStime_interval = (min([elem._SYStime for elem in self]),max([elem._SYStime for elem in self]))
            self._GPStime_interval = (min([elem._GPStime for elem in self]),max([elem._GPStime for elem in self]))

            self._lat_interval = (min([elem._lat for elem in self]),max([elem._lat for elem in self]))
            self._lon_interval = (min([elem._lon for elem in self]),max([elem._lon for elem in self]))
            self._alt_interval = (min([elem._alt for elem in self]),max([elem._alt for elem in self]))
            self._gr_spd_interval = (min([elem._gr_spd for elem in self]),max([elem._gr_spd for elem in self]))
            self._vdop_interval = (min([elem._vdop for elem in self]),max([elem._vdop for elem in self]))
            self._hdop_interval = (min([elem._hdop for elem in self]),max([elem._hdop for elem in self]))
            self._pdop_interval = (min([elem._pdop for elem in self]),max([elem._pdop for elem in self]))

            self._dopp_vel_N_interval = (min([elem._dopp_vel_N for elem in self]),max([elem._dopp_vel_N for elem in self]))
            self._dopp_vel_E_interval = (min([elem._dopp_vel_D for elem in self]),max([elem._dopp_vel_D for elem in self]))
            self._dopp_vel_D_interval = (min([elem._dopp_vel_E for elem in self]),max([elem._dopp_vel_E for elem in self]))
            self._heading_interval = (min([elem._heading for elem in self]),max([elem._heading for elem in self]))

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

        if 'gr_spd' in kwarg:
            gr_spd_i = kwarg['gr_spd'] if (len(kwarg['gr_spd']) == 2) else (kwarg['gr_spd'],kwarg['gr_spd'])
        else:
            gr_spd_i = self._gr_spd_interval

        if 'vdop' in kwarg:
            vdop_i = kwarg['vdop'] if (len(kwarg['vdop']) == 2) else (kwarg['vdop'],kwarg['vdop'])
        else:
            vdop_i = self._vdop_interval

        if 'hdop' in kwarg:
            hdop_i = kwarg['hdop'] if (len(kwarg['hdop']) == 2) else (kwarg['hdop'],kwarg['hdop'])
        else:
            hdop_i = self._hdop_interval

        if 'pdop' in kwarg:
            pdop_i = kwarg['pdop'] if (len(kwarg['pdop']) == 2) else (kwarg['pdop'],kwarg['pdop'])
        else:
            pdop_i = self._pdop_interval

        if 'dopp_vel_N' in kwarg:
            dopp_vel_N_i = kwarg['dopp_vel_N'] if (len(kwarg['dopp_vel_N']) == 2) else (kwarg['dopp_vel_N'],kwarg['dopp_vel_N'])
        else:
            dopp_vel_N_i = self._dopp_vel_N_interval

        if 'dopp_vel_E' in kwarg:
            dopp_vel_E_i = kwarg['dopp_vel_E'] if (len(kwarg['dopp_vel_E']) == 2) else (kwarg['dopp_vel_E'],kwarg['dopp_vel_E'])
        else:
            dopp_vel_E_i = self._dopp_vel_E_interval

        if 'dopp_vel_D' in kwarg:
            dopp_vel_D_i = kwarg['dopp_vel_D'] if (len(kwarg['dopp_vel_D']) == 2) else (kwarg['dopp_vel_D'],kwarg['dopp_vel_D'])
        else:
            dopp_vel_D_i = self._dopp_vel_D_interval

        if 'heading' in kwarg:
            heading_i = kwarg['heading'] if (len(kwarg['heading']) == 2) else (kwarg['heading'],kwarg['heading'])
        else:
            heading_i = self._heading_interval

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
            gr_spd_i = kwarg['selection']['gr_spd_tp'] if kwarg['selection']['gr_spd_tp'] else self._gr_spd_interval

            dopp_vel_N_i = kwarg['selection']['dopp_vel_N_tp'] if kwarg['selection']['dopp_vel_N_tp'] else self._dopp_vel_N_interval
            dopp_vel_E_i = kwarg['selection']['dopp_vel_E_tp'] if kwarg['selection']['dopp_vel_E_tp'] else self._dopp_vel_E_interval
            dopp_vel_D_i = kwarg['selection']['dopp_vel_D_tp'] if kwarg['selection']['dopp_vel_D_tp'] else self._dopp_vel_D_interval
            heading_i = kwarg['selection']['heading_tp'] if kwarg['selection']['heading_tp'] else self._heading_interval

            pdop_i = kwarg['selection']['pdop_tp'] if kwarg['selection']['pdop_tp'] else self._pdop_interval
            vdop_i = kwarg['selection']['vdop_tp'] if kwarg['selection']['vdop_tp'] else self._vdop_interval
            hdop_i = kwarg['selection']['hdop_tp'] if kwarg['selection']['hdop_tp'] else self._hdop_interval

        cnt_sel = [elem._cnt for elem in self if         (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        SYStime_sel = [elem._SYStime for elem in self if (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        GPStime_sel = [elem._GPStime for elem in self if (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        fix_val_sel = [elem._fix_val for elem in self if (  SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        n_o_satellites_sel = [elem._n_o_satellites for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        lat_sel = [elem._lat for elem in self if        (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        lon_sel = [elem._lon for elem in self if        (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        alt_sel = [elem._alt for elem in self if        (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        gr_spd_sel = [elem._gr_spd for elem in self if        (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        dopp_vel_N_sel = [elem._dopp_vel_N for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        dopp_vel_E_sel = [elem._dopp_vel_E for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        dopp_vel_D_sel = [elem._dopp_vel_D for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        heading_sel = [elem._heading for elem in self if (
                                                            SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        hdop_sel = [elem._hdop for elem in self if      (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        vdop_sel = [elem._vdop for elem in self if      (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        pdop_sel = [elem._pdop for elem in self if      (   SYStime_i[0] <= elem._SYStime <= SYStime_i[1] and
                                                            GPStime_i[0] <= elem._GPStime <= GPStime_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            fix_val_i[0] <= elem._fix_val <= fix_val_i[1] and
                                                            n_o_satellites_i[0] <= elem._n_o_satellites <= n_o_satellites_i[1] and
                                                            lat_i[0] <= elem._lat <= lat_i[1] and
                                                            lon_i[0] <= elem._lon <= lon_i[1] and
                                                            alt_i[0] <= elem._alt <= alt_i[1] and
                                                            gr_spd_i[0] <= elem._gr_spd <= gr_spd_i[1] and
                                                            dopp_vel_N_i[0] <= elem._dopp_vel_N <= dopp_vel_N_i[1] and
                                                            dopp_vel_E_i[0] <= elem._dopp_vel_E <= dopp_vel_E_i[1] and
                                                            dopp_vel_D_i[0] <= elem._dopp_vel_D <= dopp_vel_D_i[1] and
                                                            heading_i[0] <= elem._heading <= heading_i[1] and
                                                            hdop_i[0] <= elem._hdop <= hdop_i[1] and
                                                            pdop_i[0] <= elem._pdop <= pdop_i[1] and
                                                            vdop_i[0] <= elem._vdop <= vdop_i[1])]

        GPS_data = {  "SYStime": np.array(SYStime_sel),
                      "GPStime": np.array(GPStime_sel),
                      "cnt": np.array(cnt_sel),
                      "fix_val": np.array(fix_val_sel),
                      "n_o_satellites": np.array(n_o_satellites_sel),
                      "lat": np.array(lat_sel),
                      "lon": np.array(lon_sel),
                      "alt": np.array(alt_sel),
                      "gr_spd": np.array(gr_spd_sel),
                      "dopp_vel_N": np.array(dopp_vel_N_sel),
                      "dopp_vel_E": np.array(dopp_vel_E_sel),
                      "dopp_vel_D": np.array(dopp_vel_D_sel),
                      "heading": np.array(heading_sel),
                      "vdop": np.array(vdop_sel),
                      "hdop": np.array(hdop_sel),
                      "pdop": np.array(pdop_sel)}

        return GPS_data

