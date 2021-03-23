#!/usr/bin/env python
# Copyright (C) 2020-21 Andy Aschwanden

import numpy as np
from netCDF4 import Dataset as NC
from argparse import ArgumentParser


# set up the option parser
parser = ArgumentParser()
parser.add_argument("FILE", nargs="*")
parser.add_argument(
    "-s", "--scaling_factor", dest="scaling_factor", type=float, help="Scales the calving rate", default=1
)

options = parser.parse_args()
args = options.FILE
scaling_factor = options.scaling_factor

if len(args) == 0:
    nc_outfile = "seasonal_calving.nc"
elif len(args) == 1:
    nc_outfile = args[0]
else:
    print("wrong number arguments, 0 or 1 arguments accepted")
    parser.print_help()
    import sys

    sys.exit(0)


# Create netCDF file
nc = NC(nc_outfile, "w", format="NETCDF4")

nc.createDimension("time")
nc.createDimension("nb", size=2)

time = np.arange(0, 365, 1) + 0.5

h = 100

var = "time"
var_out = nc.createVariable(var, "d", dimensions=("time"))
var_out.axis = "T"
var_out.units = "days since 1980-1-1"
var_out.calendar = "365_day"
var_out.long_name = "time"
var_out.bounds = "time_bounds"
var_out[:] = time

var = "time_bounds"
var_out = nc.createVariable(var, "d", dimensions=("time", "nb"))
var_out.bounds = "time_bounds"
var_out[:, 0] = time - 0.5
var_out[:, 1] = time + 0.5


var = "frac_frac_calving_rate"
var_out = nc.createVariable(var, "f", dimensions=("time"))
var_out.units = "N m-1"

frac_calving_rate_max = 1

winter_a = 300
winter_e = 90
spring_e = 105

winter_a = 0
winter_e = 150
spring_e = 170


frac_calving_rate = np.zeros(len(time))
for k, t in enumerate(time):
    if (t < winter_e) and (t > winter_a):
        frac_calving_rate[k] = frac_calving_rate_max / np.sqrt(150) * np.sqrt(np.mod(t, 365))
    elif (t > winter_e) and (t < spring_e):
        frac_calving_rate[k] = frac_calving_rate_max - (frac_calving_rate_max / np.sqrt(20)) * np.sqrt(
            np.mod(t - winter_e, 365)
        )
    else:
        frac_calving_rate[k] = 0


for k, t in enumerate(time):
    if (t < winter_e) and (t > winter_a):
        frac_calving_rate[k] = frac_calving_rate_max - frac_calving_rate_max / np.sqrt(150) * np.sqrt(np.mod(t, 365))
    elif (t > winter_e) and (t < spring_e):
        frac_calving_rate[k] = (frac_calving_rate_max / np.sqrt(20)) * np.sqrt(np.mod(t - winter_e, 365))
    else:
        frac_calving_rate[k] = 1

var_out[:] = np.roll(frac_calving_rate, -90) * scaling_factor

nc.close()


import pylab as plt
import datetime


def set_size(w, h, ax=None):
    """ w, h: width, height in inches """

    if not ax:
        ax = plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w) / (r - l)
    figh = float(h) / (t - b)
    ax.figure.set_size_inches(figw, figh)


fontsize = 6
lw = 0.65
aspect_ratio = 0.35
markersize = 2
fig_width = 3.1  # inch
fig_height = aspect_ratio * fig_width  # inch
fig_size = [fig_width, fig_height]

params = {
    "backend": "ps",
    "axes.linewidth": 0.25,
    "lines.linewidth": lw,
    "axes.labelsize": fontsize,
    "font.size": fontsize,
    "xtick.direction": "in",
    "xtick.labelsize": fontsize,
    "xtick.major.size": 2.5,
    "xtick.major.width": 0.25,
    "ytick.direction": "in",
    "ytick.labelsize": fontsize,
    "ytick.major.size": 2.5,
    "ytick.major.width": 0.25,
    "legend.fontsize": fontsize,
    "lines.markersize": markersize,
    "font.size": fontsize,
    "figure.figsize": fig_size,
}

plt.rcParams.update(params)

positions = np.cumsum([0, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30])
labels = ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(time, frac_calving_rate)
ax.set_ylim(-0.01, 1.1)
ax.set_xlim(0, 365)
plt.xticks(positions, labels)
plt.yticks([0, frac_calving_rate_max], [0, "Max"])
ax.set_xlabel("Time [months]")
ax.set_ylabel("Calving Rate Fraction\n[1]")
set_size(3.2, 1.0)
fig.savefig("jib_seasonal_calving.pdf")
