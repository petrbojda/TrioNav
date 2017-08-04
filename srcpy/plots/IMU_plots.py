#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid.axes_size as Size
from mpl_toolkits.axes_grid import Divider
import matplotlib
import copy


def static_plot_IMU(lst_IMU, selection, fname_det,title):
    ###### Plot starts here:


    if fname_det:
        f1 = plt.figure(1, (23, 13), dpi=300)
    else:
        f1 = plt.figure(1, (15, 8))

    f1.clf()
    f1ax1 = f1.add_subplot(311)
    f1ax1.grid(True)
    plt.title(title, loc='left')
    f1ax2 = f1.add_subplot(312, sharex=f1ax1)
    f1ax2.grid(True)
    f1ax3 = f1.add_subplot(313, sharex=f1ax1)
    f1ax3.grid(True)
    # f1ax1.axis([-40, 100, -80, 80])

    #################### Left radar plot
    if lst_IMU:
        Plot_data = lst_IMU.get_array_data_sel(selection=selection)
        f1ax1.plot(Plot_data["time"], Plot_data["rotX"], label='rotation X')
        f1ax2.plot(Plot_data["time"], Plot_data["rotY"], label='rotation Y')
        f1ax3.plot(Plot_data["time"], Plot_data["rotZ"], label='rotation Z')

        plt.draw()

    f1ax3.set_xlabel('time [ms]')
    f1ax1.set_ylabel('rate X [deg/sec]')
    f1ax2.set_ylabel('rate Y [deg/sec]')
    f1ax3.set_ylabel('rate Z [deg/sec]')

    # tit = "In Selected Beams: L=%d R=%d" % (number_of_dets_left_processed, number_of_dets_right_processed)
    # f1.suptitle(tit, fontsize=14, fontweight='bold')

    if fname_det:
        f1.savefig(fname_det)
    else:
        plt.show()

