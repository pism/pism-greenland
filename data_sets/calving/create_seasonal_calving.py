#!/usr/bin/env python
# Copyright (C) 2020-21 Andy Aschwanden

from argparse import ArgumentParser
from calendar import isleap
from dateutil import rrule
from datetime import datetime
import numpy as np
from netCDF4 import Dataset as NC
import cftime

# set up the option parser
parser = ArgumentParser()
parser.add_argument("FILE", nargs="*")
parser.add_argument(
    "-s",
    "--scaling_factor",
    dest="scaling_factor",
    type=float,
    help="Scales the calving rate",
    default=1,
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

start_year = 1980
end_year = 2021
start_date = datetime(start_year, 1, 1)
end_date = datetime(end_year, 1, 1)

calendar = "standard"
units = "days since 1980-1-1"

sampling_interval = "daily"
rd = {
    "daily": rrule.DAILY,
    "weekly": rrule.WEEKLY,
    "monthly": rrule.MONTHLY,
    "yearly": rrule.YEARLY,
}

bnds_datelist = list(
    rrule.rrule(rd[sampling_interval], dtstart=start_date, until=end_date)
)
# calculate the days since refdate, including refdate, with time being the
bnds_interval_since_refdate = cftime.date2num(bnds_datelist, units, calendar=calendar)
time_interval_since_refdate = (
    bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2
)

# mid-point value:
# time[n] = (bnds[n] + bnds[n+1]) / 2
time_interval_since_refdate = (
    bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2
)

nt = len(time_interval_since_refdate)

# Create netCDF file
nc = NC(nc_outfile, "w", format="NETCDF4", compression_level=2)

nc.createDimension("time")
nc.createDimension("nb", size=2)

var = "time"
var_out = nc.createVariable(var, "d", dimensions=("time"))
var_out.axis = "T"
var_out.units = "days since 1980-1-1"
var_out.long_name = "time"
var_out.bounds = "time_bounds"
var_out[:] = time_interval_since_refdate

var = "time_bounds"
var_out = nc.createVariable(var, "d", dimensions=("time", "nb"))
var_out.bounds = "time_bounds"
var_out[:, 0] = bnds_interval_since_refdate[0:-1]
var_out[:, 1] = bnds_interval_since_refdate[1::]


var = "frac_calving_rate"
var_out = nc.createVariable(var, "f", dimensions=("time"))
var_out.units = "1"

frac_calving_rate_max = 1

winter_a = 300
winter_e = 90
spring_e = 105

winter_a = 0
winter_e = 150
spring_e = 170

idx = 0
for year in range(start_year, end_year):
    print(f"Preparing Year {year}")
    if isleap(year):
        year_length = 366
    else:
        year_length = 365

    frac_calving_rate = np.zeros(year_length)
    for t in range(year_length):
        if (t <= winter_e) and (t >= winter_a):
            frac_calving_rate[
                t
            ] = frac_calving_rate_max - frac_calving_rate_max / np.sqrt(
                winter_e
            ) * np.sqrt(
                np.mod(t, year_length)
            )
        elif (t > winter_e) and (t < spring_e):
            frac_calving_rate[t] = (
                frac_calving_rate_max / np.sqrt(spring_e - winter_e)
            ) * np.sqrt(np.mod(t - winter_e, year_length))
        else:
            frac_calving_rate[t] = 1
            if year > 2001:
                frac_calving_rate[t] = (
                    frac_calving_rate_max / np.sqrt(spring_e - winter_e)
                ) * np.sqrt(np.mod(t - winter_e, year_length))

    frac_calving_rate = np.roll(frac_calving_rate, -90) * scaling_factor
    var_out[idx::] = frac_calving_rate
    idx += year_length

nc.close()


import pylab as plt
import datetime


def set_size(w, h, ax=None):
    """w, h: width, height in inches"""

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
labels = [
    "Oct",
    "Nov",
    "Dec",
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(range(len(frac_calving_rate)), np.roll(frac_calving_rate, 90))
ax.set_ylim(-0.01, 3.1)
ax.set_xlim(0, len(frac_calving_rate))
plt.xticks(positions, labels)
plt.yticks([0, frac_calving_rate_max], [0, "Max"])
ax.set_xlabel("Time [months]")
ax.set_ylabel("Calving Rate Fraction\n[1]")
set_size(3.2, 1.0)
fig.savefig(f"jib_seasonal_calving_{scaling_factor}_{start_year}_{end_year}.pdf")
