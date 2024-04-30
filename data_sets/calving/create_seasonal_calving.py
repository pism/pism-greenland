#!/usr/bin/env python
# Copyright (C) 2020-23 Andy Aschwanden

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
    "--year_high",
    type=float,
    help="Start when high values are applied.",
    default=2001,
)
parser.add_argument(
    "--calving_low",
    type=float,
    help="Start when high values are applied.",
    default=1.0,
)
parser.add_argument(
    "--calving_high",
    type=float,
    help="Start when high values are applied.",
    default=1.5,
)

options = parser.parse_args()
args = options.FILE
r_low = options.calving_low
r_high = options.calving_high
year_high = options.year_high

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


winter_a = 300
winter_e = 90
spring_e = 105

winter_a = 0
winter_e = 150
spring_e = 170

def annual_calving(year_length, frac_calving_rate_max):
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
            frac_calving_rate[t] = frac_calving_rate_max
    return frac_calving_rate

idx = 0
for year in range(start_year, end_year):
    print(f"Preparing Year {year}")
    if isleap(year):
        year_length = 366
    else:
        year_length = 365

    if year <= year_high:
        frac_calving_rate = annual_calving(year_length, r_low)
    else:
        frac_calving_rate = annual_calving(year_length, r_high)

    frac_calving_rate = np.roll(frac_calving_rate, -90) 
    var_out[idx::] = frac_calving_rate
    idx += year_length

nc.close()

