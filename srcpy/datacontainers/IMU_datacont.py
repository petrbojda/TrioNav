import scipy.io as sio
import numpy as np




class MeasuredPoint_IMU(object):
    def __init__(self, time=None, cnt=None, IMU_accX=None, IMU_DaccX=None, IMU_accY=None, IMU_DaccY=None, IMU_accZ=None, IMU_DaccZ=None,
                 IMU_rotX=None, IMU_DrotX=None, IMU_rotY=None, IMU_DrotY=None, IMU_rotZ=None, IMU_DrotZ=None):

        self._time = time
        self._cnt = cnt
        self._IMU_accX = IMU_accX
        self._IMU_DaccX =IMU_DaccX
        self._IMU_accY = IMU_accY
        self._IMU_DaccY =IMU_DaccY
        self._IMU_accZ = IMU_accZ
        self._IMU_DaccZ =IMU_DaccZ

        self._IMU_rotX = IMU_rotX
        self._IMU_DrotX =IMU_DrotX
        self._IMU_rotY = IMU_rotY
        self._IMU_DrotY =IMU_DrotY
        self._IMU_rotZ = IMU_rotZ
        self._IMU_DrotZ =IMU_DrotZ


class List_MP_IMU(list):

    def __init__(self,IMU_sampling_rate = 200):
        super().__init__()
        self._sampling_rate = IMU_sampling_rate
        self._time_interval = (0,0)

        self._IMU_accX_interval = (0,0)
        self._IMU_accY_interval = (0,0)
        self._IMU_accZ_interval = (0,0)
        self._IMU_DaccX_interval = (0,0)
        self._IMU_DaccY_interval = (0,0)
        self._IMU_DaccZ_interval = (0,0)

        self._IMU_rotX_interval = (0,0)
        self._IMU_rotY_interval = (0,0)
        self._IMU_rotZ_interval = (0,0)
        self._IMU_DrotX_interval = (0,0)
        self._IMU_DrotY_interval = (0,0)
        self._IMU_DrotZ_interval = (0,0)


    def append_from_m_file(self, data_path):
        measured_data_mfile = sio.loadmat(data_path)
        raw_data = measured_data_mfile["data"]
        no_d = len(raw_data)
        for itr in range(0, no_d ):
            self.append(MeasuredPoint_IMU ( cnt = int(raw_data[itr, 0]),
                                            time = int(raw_data[itr, 0]) / self._sampling_rate,
                                            IMU_accX=float(raw_data[itr, 2]),
                                            IMU_accY=float(raw_data[itr, 4]),
                                            IMU_accZ=float(raw_data[itr, 6]),
                                            IMU_DaccX=float(raw_data[itr, 3]),
                                            IMU_DaccY=float(raw_data[itr, 5]),
                                            IMU_DaccZ=float(raw_data[itr, 7]),
                                            IMU_rotX=float(raw_data[itr, 10]),
                                            IMU_rotY=float(raw_data[itr, 12]),
                                            IMU_rotZ=float(raw_data[itr, 14]),
                                            IMU_DrotX=float(raw_data[itr, 11]),
                                            IMU_DrotY=float(raw_data[itr, 13]),
                                            IMU_DrotZ=float(raw_data[itr, 15])))


        self._cnt_interval = (min([elem._cnt for elem in self]),max([elem._cnt for elem in self]))
        self._time_interval = (min([elem._time for elem in self]),max([elem._time for elem in self]))

        self._IMU_accX_interval = (min([elem._IMU_accX for elem in self]),max([elem._IMU_accX for elem in self]))
        self._IMU_accY_interval = (min([elem._IMU_accY for elem in self]),max([elem._IMU_accY for elem in self]))
        self._IMU_accZ_interval = (min([elem._IMU_accZ for elem in self]),max([elem._IMU_accZ for elem in self]))
        self._IMU_DaccX_interval = (min([elem._IMU_DaccX for elem in self]),max([elem._IMU_DaccX for elem in self]))
        self._IMU_DaccY_interval = (min([elem._IMU_DaccY for elem in self]),max([elem._IMU_DaccY for elem in self]))
        self._IMU_DaccZ_interval = (min([elem._IMU_DaccZ for elem in self]),max([elem._IMU_DaccZ for elem in self]))

        self._IMU_rotX_interval = (min([elem._IMU_rotX for elem in self]),max([elem._IMU_rotX for elem in self]))
        self._IMU_rotY_interval = (min([elem._IMU_rotY for elem in self]),max([elem._IMU_rotY for elem in self]))
        self._IMU_rotZ_interval = (min([elem._IMU_rotZ for elem in self]),max([elem._IMU_rotZ for elem in self]))
        self._IMU_DrotX_interval = (min([elem._IMU_DrotX for elem in self]),max([elem._IMU_DrotX for elem in self]))
        self._IMU_DrotY_interval = (min([elem._IMU_DrotY for elem in self]),max([elem._IMU_DrotY for elem in self]))
        self._IMU_DrotZ_interval = (min([elem._IMU_DrotZ for elem in self]),max([elem._IMU_DrotZ for elem in self]))

    def append_from_dict(self, dictIMU):

        no_d = len(dictIMU["cnt"])
        for itr in range(0, no_d ):
            self.append(MeasuredPoint_IMU ( cnt = dictIMU["cnt"][itr],
                                            time = dictIMU["time"][itr],
                                            IMU_accX = dictIMU["accX"][itr],
                                            IMU_accY=dictIMU["accY"][itr],
                                            IMU_accZ=dictIMU["accZ"][itr],
                                            IMU_DaccX=dictIMU["DaccX"][itr],
                                            IMU_DaccY=dictIMU["DaccY"][itr],
                                            IMU_DaccZ=dictIMU["DaccZ"][itr],
                                            IMU_rotX=dictIMU["rotX"][itr],
                                            IMU_rotY=dictIMU["rotY"][itr],
                                            IMU_rotZ=dictIMU["rotZ"][itr],
                                            IMU_DrotX=dictIMU["DrotX"][itr],
                                            IMU_DrotY=dictIMU["DrotY"][itr],
                                            IMU_DrotZ=dictIMU["DrotZ"][itr]))


        self._cnt_interval = (min([elem._cnt for elem in self]),max([elem._cnt for elem in self]))
        self._time_interval = (min([elem._time for elem in self]),max([elem._time for elem in self]))

        self._IMU_accX_interval = (min([elem._IMU_accX for elem in self]),max([elem._IMU_accX for elem in self]))
        self._IMU_accY_interval = (min([elem._IMU_accY for elem in self]),max([elem._IMU_accY for elem in self]))
        self._IMU_accZ_interval = (min([elem._IMU_accZ for elem in self]),max([elem._IMU_accZ for elem in self]))
        self._IMU_DaccX_interval = (min([elem._IMU_DaccX for elem in self]),max([elem._IMU_DaccX for elem in self]))
        self._IMU_DaccY_interval = (min([elem._IMU_DaccY for elem in self]),max([elem._IMU_DaccY for elem in self]))
        self._IMU_DaccZ_interval = (min([elem._IMU_DaccZ for elem in self]),max([elem._IMU_DaccZ for elem in self]))

        self._IMU_rotX_interval = (min([elem._IMU_rotX for elem in self]),max([elem._IMU_rotX for elem in self]))
        self._IMU_rotY_interval = (min([elem._IMU_rotY for elem in self]),max([elem._IMU_rotY for elem in self]))
        self._IMU_rotZ_interval = (min([elem._IMU_rotZ for elem in self]),max([elem._IMU_rotZ for elem in self]))
        self._IMU_DrotX_interval = (min([elem._IMU_DrotX for elem in self]),max([elem._IMU_DrotX for elem in self]))
        self._IMU_DrotY_interval = (min([elem._IMU_DrotY for elem in self]),max([elem._IMU_DrotY for elem in self]))
        self._IMU_DrotZ_interval = (min([elem._IMU_DrotZ for elem in self]),max([elem._IMU_DrotZ for elem in self]))


    def get_count_interval(self):
        self._cnt_interval = (min([elem._cnt for elem in self]),max([elem._cnt for elem in self]))
        return self._cnt_interval

    def get_time_interval(self):
        self._time_interval = (min([elem._time for elem in self]),max([elem._time for elem in self]))
        return self._time_interval

    def get_time_duration(self):
        time_duration = (self._time_interval[1] - self._time_interval[0])
        return time_duration

    def get_IMU_data_attr(self):
        time_duration = (self._time_interval[1] - self._time_interval[0])
        return time_duration

    def get_array_data_sel(self, **kwarg):

        if 'cnt' in kwarg:
            cnt_i = kwarg['cnt'] if (len(kwarg['cnt']) == 2) else (kwarg['cnt'],kwarg['cnt'])
        else:
            cnt_i = self._cnt_interval

        if 'time' in kwarg:
            time_i = kwarg['time'] if (len(kwarg['time']) == 2) else (kwarg['time'],kwarg['time'])
        else:
            time_i = self._time_interval

        if 'rotX' in kwarg:
            rotX_i = kwarg['rotX'] if (len(kwarg['rotX']) == 2) else (kwarg['rotX'],kwarg['rotX'])
        else:
            rotX_i = self._IMU_rotX_interval

        if 'rotY' in kwarg:
            rotY_i = kwarg['rotY'] if (len(kwarg['rotY']) == 2) else (kwarg['rotY'],kwarg['rotY'])
        else:
            rotY_i = self._IMU_rotY_interval

        if 'rotZ' in kwarg:
            rotZ_i = kwarg['rotZ'] if (len(kwarg['rotZ']) == 2) else (kwarg['rotZ'],kwarg['rotZ'])
        else:
            rotZ_i = self._IMU_rotZ_interval

        if 'DrotX' in kwarg:
            DrotX_i = kwarg['DrotX'] if (len(kwarg['DrotX']) == 2) else (kwarg['DrotX'],kwarg['DrotX'])
        else:
            DrotX_i = self._IMU_DrotX_interval

        if 'DrotY' in kwarg:
            DrotY_i = kwarg['DrotY'] if (len(kwarg['DrotY']) == 2) else (kwarg['DrotY'],kwarg['DrotY'])
        else:
            DrotY_i = self._IMU_DrotY_interval

        if 'DrotZ' in kwarg:
            DrotZ_i = kwarg['DrotZ'] if (len(kwarg['DrotZ']) == 2) else (kwarg['DrotZ'],kwarg['DrotZ'])
        else:
            DrotZ_i = self._IMU_DrotZ_interval

        if 'accX' in kwarg:
            accX_i = kwarg['accX'] if (len(kwarg['accX']) == 2) else (kwarg['accX'],kwarg['accX'])
        else:
            accX_i = self._IMU_accX_interval

        if 'accY' in kwarg:
            accY_i = kwarg['accY'] if (len(kwarg['accY']) == 2) else (kwarg['accY'],kwarg['accY'])
        else:
            accY_i = self._IMU_accY_interval

        if 'accZ' in kwarg:
            accZ_i = kwarg['accZ'] if (len(kwarg['accZ']) == 2) else (kwarg['accZ'],kwarg['accZ'])
        else:
            accZ_i = self._IMU_accZ_interval

        if 'DaccX' in kwarg:
            DaccX_i = kwarg['DaccX'] if (len(kwarg['DaccX']) == 2) else (kwarg['DaccX'],kwarg['DaccX'])
        else:
            DaccX_i = self._IMU_DaccX_interval

        if 'DaccY' in kwarg:
            DaccY_i = kwarg['DaccY'] if (len(kwarg['DaccY']) == 2) else (kwarg['DaccY'],kwarg['DaccY'])
        else:
            DaccY_i = self._IMU_DaccY_interval

        if 'DaccZ' in kwarg:
            DaccZ_i = kwarg['DaccZ'] if (len(kwarg['DaccZ']) == 2) else (kwarg['DaccZ'],kwarg['DaccZ'])
        else:
            DaccZ_i = self._IMU_DaccZ_interval


        if 'selection' in kwarg:
            rotX_i = kwarg['selection']['rotX_tp'] if kwarg['selection']['rotX_tp'] else self._IMU_rotX_interval
            rotY_i = kwarg['selection']['rotY_tp'] if kwarg['selection']['rotY_tp'] else self._IMU_rotY_interval
            rotZ_i = kwarg['selection']['rotZ_tp'] if kwarg['selection']['rotZ_tp'] else self._IMU_rotZ_interval

            DrotX_i = kwarg['selection']['DrotX_tp'] if kwarg['selection']['DrotX_tp'] else self._IMU_DrotX_interval
            DrotY_i = kwarg['selection']['DrotY_tp'] if kwarg['selection']['DrotY_tp'] else self._IMU_DrotY_interval
            DrotZ_i = kwarg['selection']['DrotZ_tp'] if kwarg['selection']['DrotZ_tp'] else self._IMU_DrotZ_interval

            accX_i = kwarg['selection']['accX_tp'] if kwarg['selection']['accX_tp'] else self._IMU_accX_interval
            accY_i = kwarg['selection']['accY_tp'] if kwarg['selection']['accY_tp'] else self._IMU_accY_interval
            accZ_i = kwarg['selection']['accZ_tp'] if kwarg['selection']['accZ_tp'] else self._IMU_accZ_interval

            DaccX_i = kwarg['selection']['DaccX_tp'] if kwarg['selection']['DaccX_tp'] else self._IMU_DaccX_interval
            DaccY_i = kwarg['selection']['DaccY_tp'] if kwarg['selection']['DaccY_tp'] else self._IMU_DaccY_interval
            DaccZ_i = kwarg['selection']['DaccZ_tp'] if kwarg['selection']['DaccZ_tp'] else self._IMU_DaccZ_interval

            cnt_i = kwarg['selection']['cnt_tp'] if kwarg['selection']['cnt_tp'] else self._cnt_interval

            time_i = kwarg['selection']['time_tp'] if kwarg['selection']['time_tp'] else self._time_interval


        rotX_sel = [elem._IMU_rotX for elem in self if (    time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        rotY_sel = [elem._IMU_rotY for elem in self if (    time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        rotZ_sel = [elem._IMU_rotZ for elem in self if (    time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        DrotX_sel = [elem._IMU_DrotX for elem in self if (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        DrotY_sel = [elem._IMU_DrotY for elem in self if (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        DrotZ_sel = [elem._IMU_DrotZ for elem in self if (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        accX_sel = [elem._IMU_accX for elem in self if (    time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        accY_sel = [elem._IMU_accY for elem in self if (    time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        accZ_sel = [elem._IMU_accZ for elem in self if (    time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        DaccX_sel = [elem._IMU_DaccX for elem in self if (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        DaccY_sel = [elem._IMU_DaccY for elem in self if (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        DaccZ_sel = [elem._IMU_DaccZ for elem in self if (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        cnt_sel = [elem._cnt for elem in self if         (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        time_sel = [elem._time for elem in self if       (  time_i[0] <= elem._time <= time_i[1] and
                                                            cnt_i[0] <= elem._cnt <= cnt_i[1] and
                                                            rotX_i[0] <= elem._IMU_rotX <= rotX_i[1] and
                                                            rotY_i[0] <= elem._IMU_rotY <= rotY_i[1] and
                                                            rotZ_i[0] <= elem._IMU_rotZ <= rotZ_i[1] and
                                                            DrotX_i[0] <= elem._IMU_DrotX <= DrotX_i[1] and
                                                            DrotY_i[0] <= elem._IMU_DrotY <= DrotY_i[1] and
                                                            DrotZ_i[0] <= elem._IMU_DrotZ <= DrotZ_i[1] and
                                                            accX_i[0] <= elem._IMU_accX <= accX_i[1] and
                                                            accY_i[0] <= elem._IMU_accY <= accY_i[1] and
                                                            accZ_i[0] <= elem._IMU_accZ <= accZ_i[1] and
                                                            DaccX_i[0] <= elem._IMU_DaccX <= DaccX_i[1] and
                                                            DaccY_i[0] <= elem._IMU_DaccY <= DaccY_i[1] and
                                                            DaccZ_i[0] <= elem._IMU_DaccZ <= DaccZ_i[1])]

        IMU_data = {  "rotX": np.array(rotX_sel),
                      "rotY": np.array(rotY_sel),
                      "rotZ": np.array(rotZ_sel),
                      "DrotX": np.array(DrotX_sel),
                      "DrotY": np.array(DrotY_sel),
                      "DrotZ": np.array(DrotZ_sel),
                      "accX": np.array(accX_sel),
                      "accY": np.array(accY_sel),
                      "accZ": np.array(accZ_sel),
                      "DaccX": np.array(DaccX_sel),
                      "DaccY": np.array(DaccY_sel),
                      "DaccZ": np.array(DaccZ_sel),
                      "time": np.array(time_sel),
                      "cnt": np.array(cnt_sel)}

        return IMU_data



