#!/usr/bin/env python

#import numpy as np
import matplotlib.pyplot as plt
#import mpl_toolkits.axes_grid.axes_size as Size
#from mpl_toolkits.axes_grid import Divider
#import matplotlib
#import copy


def all_raw(lst_RTK, selection, fname, title, fig):

    if fname:
        f1 = plt.figure(fig + 0, (23, 13), dpi=300)
        f2 = plt.figure(fig + 1, (23, 13), dpi=300)
        f3 = plt.figure(fig + 2, (23, 13), dpi=300)
        f4 = plt.figure(fig + 3, (23, 13), dpi=300)
    else:
        f1 = plt.figure(fig + 0, (15, 8))
        f2 = plt.figure(fig + 1, (15, 8))
        f3 = plt.figure(fig + 2, (15, 8))
        f4 = plt.figure(fig + 3, (15, 8))

    tit_ecef = 'ECEF ' + title
    tit_stdNEU = 'STD of North, East and Up ' + title
    tit_RTKage = 'RTK Age and Ratio ' + title
    tit_stdNeEuUn = 'STD of N-E, E-U and U-N ' + title

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

    # std N, E, U plotted to f2
    f2.clf()
    f2ax1 = f2.add_subplot(311, sharex=f1ax1)
    f2ax1.grid(True)
    #plt.title(tit_stdNEU, loc='left')
    f2.suptitle(tit_stdNEU, fontsize=14, fontweight='bold')
    f2ax2 = f2.add_subplot(312, sharex=f1ax1)
    f2ax2.grid(True)
    f2ax3 = f2.add_subplot(313, sharex=f1ax1)
    f2ax3.grid(True)

    # RTK age + RTK ratio plotted to f3
    f3.clf()
    f3ax1 = f3.add_subplot(211, sharex=f1ax1)
    f3ax1.grid(True)
    #plt.title(tit_RTKage, loc='left')
    f3.suptitle(tit_RTKage, fontsize=14, fontweight='bold')
    f3ax2 = f3.add_subplot(212, sharex=f1ax1)
    f3ax2.grid(True)

    # std NE, EU, UN plotted to f4
    f4.clf()
    f4ax1 = f4.add_subplot(311, sharex=f1ax1)
    f4ax1.grid(True)
    #plt.title(tit_stdNeEuUn, loc='left')
    f4.suptitle(tit_stdNeEuUn, fontsize=14, fontweight='bold')
    f4ax2 = f4.add_subplot(312, sharex=f1ax1)
    f4ax2.grid(True)
    f4ax3 = f4.add_subplot(313, sharex=f1ax1)
    f4ax3.grid(True)

    if lst_RTK:
        Plot_data = lst_RTK.get_array_data_sel(selection=selection)
        f1ax1.plot(Plot_data["SYStime"], Plot_data["lat"],'-b', label='Latitude')
        f1ax2.plot(Plot_data["SYStime"], Plot_data["lon"],'-b', label='Longitude')
        f1ax3.plot(Plot_data["SYStime"], Plot_data["alt"],'-b', label='Altitude')
        plt.draw()

        f2ax1.plot(Plot_data["SYStime"], Plot_data["stdNorth"],'-b', label='North')
        f2ax2.plot(Plot_data["SYStime"], Plot_data["stdEast"],'-b', label='East')
        f2ax3.plot(Plot_data["SYStime"], Plot_data["stdUp"],'-b', label='Down')
        plt.draw()

        f3ax1.plot(Plot_data["SYStime"], Plot_data["RTKage"],'-b', label='age')
        f3ax2.plot(Plot_data["SYStime"], Plot_data["RTKratio"],'-b', label='ratio')
        plt.draw()

        f4ax1.plot(Plot_data["SYStime"], Plot_data["stdNE"],'-b', label='stdNE')
        f4ax2.plot(Plot_data["SYStime"], Plot_data["stdEU"],'-b', label='stdEU')
        f4ax3.plot(Plot_data["SYStime"], Plot_data["stdUN"],'-b', label='stdUN')
        plt.draw()

    f1ax3.set_xlabel('time [ms]')
    f1ax1.set_ylabel('Lat [deg]')
    f1ax2.set_ylabel('Lon [deg]')
    f1ax3.set_ylabel('Alt [meters]')

    f2ax3.set_xlabel('time [ms]')
    f2ax1.set_ylabel('N [meters]')
    f2ax2.set_ylabel('E [meters]')
    f2ax3.set_ylabel('D [meters]')

    f3ax2.set_xlabel('time [ms]')
    f3ax1.set_ylabel('age [sec]')
    f3ax2.set_ylabel('ratio [-]')

    f4ax3.set_xlabel('time [ms]')
    f4ax1.set_ylabel('NE [meters]')
    f4ax2.set_ylabel('EU [meters]')
    f4ax3.set_ylabel('UN [meters]')

    if fname:
        fname_ecef = fname + '_ecef'
        fname_stdNEU = fname + '_stdNEU'
        fname_RTKage = fname + '_RTKage'
        fname_stdNeEuUn = fname + '_stdNeEuUn'
        f1.savefig(fname_ecef)
        f2.savefig(fname_stdNEU)
        f3.savefig(fname_RTKage)
        f4.savefig(fname_stdNeEuUn)
    else:
        plt.show()
