#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from netCDF4 import Dataset as NC
import numpy as np
import pylab as plt


def set_size(w, h, ax=None):

    """ 
    w, h: width, height in inches
    """

    if not ax:
        ax = plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w) / (r - l)
    figh = float(h) / (t - b)
    ax.figure.set_size_inches(figw, figh)


nc = NC("fldsum_synth_jib_runoff_g1000m.nc", "r")
water_input_rate = np.squeeze(nc.variables["water_input_rate"][:])
nc.close()

nc_r = NC("2019_12_routing/spatial/fldmax_ex_synth_jib_g1000m_id_ROUTING_0_10.nc", "r")
nc_s = NC("2019_12_steady/spatial/fldmax_ex_synth_jib_g1000m_id_STEADY_0_10.nc", "r")

calving_r = np.squeeze(nc_r.variables["vonmises_calving_rate"][:])
fm_r = np.squeeze(nc_r.variables["frontal_melt_rate"][:])
fm_s = np.squeeze(nc_s.variables["frontal_melt_rate"][:])

time = range(0, 365)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(time, water_input_rate)
ax.set_xlabel("days of year")
ax.set_ylabel("water input rate (Gt/yr")

fig.savefig("water_input_rate.pdf", bbox_inches="tight")


fig, ax = plt.subplots(2, 1, sharex="col", figsize=[4.75, 3.5])
fig.subplots_adjust(hspace=0.25, wspace=0.05)

ax[0].plot(time, water_input_rate)
ax[1].plot(time, fm_r, label="ROUTING")
ax[1].plot(time, fm_s, label="STEADY")
ax[-1].set_xlabel("days of year")
ax[0].set_ylabel("water input\n(Gt/yr)")
ax[1].set_ylabel("frontal melt rate\n(m/day)")
legend = ax[1].legend(ncol=1)
legend.get_frame().set_linewidth(0.0)
legend.get_frame().set_alpha(0.0)

fig.savefig("routing_vs_steady.pdf", bbox_inches="tight")
