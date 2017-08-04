#!/usr/bin/env python

#import numpy as np
import matplotlib.pyplot as plt
#import mpl_toolkits.axes_grid.axes_size as Size
#from mpl_toolkits.axes_grid import Divider
#import matplotlib
#import copy


def three_rotacc_raw(lst_IMU, selection, fname,title,fig):

    if fname:
        f1 = plt.figure(fig, (23, 13), dpi=300)
        f2 = plt.figure(fig + 1, (23, 13), dpi=300)
    else:
        f1 = plt.figure(fig, (15, 8))
        f2 = plt.figure(fig + 1, (15, 8))

    tit_rot = 'Angular rates ' + title
    tit_acc = 'Accelerations ' + title

    # Angular rates plotted to f1
    f1.clf()
    f1ax1 = f1.add_subplot(311)
    f1ax1.grid(True)
    #plt.title(tit_rot, loc='left')
    f1.suptitle(tit_rot, fontsize=14, fontweight='bold')
    f1ax2 = f1.add_subplot(312, sharex=f1ax1)
    f1ax2.grid(True)
    f1ax3 = f1.add_subplot(313, sharex=f1ax1)
    f1ax3.grid(True)
    # f1ax1.axis([-40, 100, -80, 80])

    # Accelerations plotted to f2
    f2.clf()
    f2ax1 = f2.add_subplot(311)
    f2ax1.grid(True)
    #plt.title(title, loc='left')
    f2.suptitle(tit_acc, fontsize=14, fontweight='bold')
    f2ax2 = f2.add_subplot(312, sharex=f2ax1)
    f2ax2.grid(True)
    f2ax3 = f2.add_subplot(313, sharex=f2ax1)
    f2ax3.grid(True)

    if lst_IMU:
        Plot_data = lst_IMU.get_array_data_sel(selection=selection)
        f1ax1.plot(Plot_data["time"], Plot_data["rotX"],'-b', label='ars X')
        f1ax2.plot(Plot_data["time"], Plot_data["rotY"],'-b', label='ars Y')
        f1ax3.plot(Plot_data["time"], Plot_data["rotZ"],'-b', label='ars Z')
        f1ax1.plot(Plot_data["time"], Plot_data["DrotX"],'-r', label='delta ars X')
        f1ax2.plot(Plot_data["time"], Plot_data["DrotY"],'-r', label='delta ars Y')
        f1ax3.plot(Plot_data["time"], Plot_data["DrotZ"],'-r', label='delta ars Z')
        plt.draw()

        f2ax1.plot(Plot_data["time"], Plot_data["accX"],'-b', label='acc X')
        f2ax2.plot(Plot_data["time"], Plot_data["accY"],'-b', label='acc Y')
        f2ax3.plot(Plot_data["time"], Plot_data["accZ"],'-b', label='acc Z')
        f2ax1.plot(Plot_data["time"], Plot_data["DaccX"],'-r', label='delta acc X')
        f2ax2.plot(Plot_data["time"], Plot_data["DaccY"],'-r', label='delta acc Y')
        f2ax3.plot(Plot_data["time"], Plot_data["DaccZ"],'-r', label='delta acc Z')
        plt.draw()

    f1ax3.set_xlabel('time [ms]')
    f1ax1.set_ylabel('rate X [deg/sec]')
    f1ax2.set_ylabel('rate Y [deg/sec]')
    f1ax3.set_ylabel('rate Z [deg/sec]')

    f2ax3.set_xlabel('time [ms]')
    f2ax1.set_ylabel('acc X [g]')
    f2ax2.set_ylabel('acc Y [g]')
    f2ax3.set_ylabel('acc Z [g]')

    if fname:
        fname_rot = fname + '_rot'
        fname_acc = fname + '_acc'
        f1.savefig(fname_rot)
        f2.savefig(fname_acc)
    else:
        plt.show()

