#!/usr/bin/env python

#import numpy as np
import matplotlib.pyplot as plt
#import mpl_toolkits.axes_grid.axes_size as Size
#from mpl_toolkits.axes_grid import Divider
#import matplotlib
#import copy


def all_raw(lst_GPS, selection, fname, title, fig):

    if fname:
        f1 = plt.figure(fig + 0, (23, 13), dpi=300)
        f2 = plt.figure(fig + 1, (23, 13), dpi=300)
        f3 = plt.figure(fig + 2, (23, 13), dpi=300)
    else:
        f1 = plt.figure(fig + 0, (15, 8))
        f2 = plt.figure(fig + 1, (15, 8))
        f3 = plt.figure(fig + 2, (15, 8))

    tit_ecef = 'ECEF ' + title
    tit_NEDvel = 'NED Velocities ' + title
    tit_head_grspd = 'Ground Projection ' + title

    # ECEF Lat, Long, Alt plotted to f1
    f1.clf()
    f1ax1 = f1.add_subplot(311)
    f1ax1.grid(True)
    #plt.title(tit_ecef, loc='left')
    f1.suptitle(tit_ecef, fontsize=14, fontweight='bold')
    f1ax2 = f1.add_subplot(312, sharex=f1ax1)
    f1ax2.grid(True)
    f1ax3 = f1.add_subplot(313, sharex=f1ax1)
    f1ax3.grid(True)
    # f1ax1.axis([-40, 100, -80, 80])

    # NED velocities plotted to f2
    f2.clf()
    f2ax1 = f2.add_subplot(311)
    f2ax1.grid(True)
    #plt.title(tit_NEDvel, loc='left')
    f2.suptitle(tit_NEDvel, fontsize=14, fontweight='bold')
    f2ax2 = f2.add_subplot(312, sharex=f2ax1)
    f2ax2.grid(True)
    f2ax3 = f2.add_subplot(313, sharex=f2ax1)
    f2ax3.grid(True)

    # Ground Speed and Heading plotted to f3
    f3.clf()
    f3ax1 = f3.add_subplot(211)
    f3ax1.grid(True)
    #plt.title(tit_head_grspd, loc='left')
    f3.suptitle(tit_head_grspd, fontsize=14, fontweight='bold')
    f3ax2 = f3.add_subplot(212, sharex=f3ax1)
    f3ax2.grid(True)

    if lst_GPS:
        Plot_data = lst_GPS.get_array_data_sel(selection=selection)
        f1ax1.plot(Plot_data["SYStime"], Plot_data["lat"],'-b', label='Latitude')
        f1ax2.plot(Plot_data["SYStime"], Plot_data["lon"],'-b', label='Longitude')
        f1ax3.plot(Plot_data["SYStime"], Plot_data["alt"],'-b', label='Altitude')
        plt.draw()

        f2ax1.plot(Plot_data["SYStime"], Plot_data["dopp_vel_N"],'-b', label='North')
        f2ax2.plot(Plot_data["SYStime"], Plot_data["dopp_vel_E"],'-b', label='East')
        f2ax3.plot(Plot_data["SYStime"], Plot_data["dopp_vel_D"],'-b', label='Down')
        plt.draw()

        f3ax1.plot(Plot_data["SYStime"], Plot_data["heading"],'-b', label='heading')
        f3ax2.plot(Plot_data["SYStime"], Plot_data["gr_spd"],'-b', label='ground speed')
        plt.draw()

    f1ax3.set_xlabel('time [ms]')
    f1ax1.set_ylabel('Lat [deg]')
    f1ax2.set_ylabel('Lon [deg]')
    f1ax3.set_ylabel('Alt [meters]')

    f2ax3.set_xlabel('time [ms]')
    f2ax1.set_ylabel('N vel [m/sec]')
    f2ax2.set_ylabel('E vel [m/sec]')
    f2ax3.set_ylabel('D vel [m/sec]')

    f3ax2.set_xlabel('time [ms]')
    f3ax1.set_ylabel('Head [deg]')
    f3ax2.set_ylabel('Gr. spd [m/sec]')

    if fname:
        fname_ecef = fname + '_ecef'
        fname_NEDvel = fname + '_NEDvel'
        fname_head_grspd = fname + '_head_grspd'
        f1.savefig(fname_ecef)
        f2.savefig(fname_NEDvel)
        f3.savefig(fname_head_grspd)
    else:
        plt.show()
