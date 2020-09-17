#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

import numpy as np
from netCDF4 import Dataset as NC
from argparse import ArgumentParser


# set up the option parser
parser = ArgumentParser()
parser.add_argument("FILE", nargs="*")
parser.add_argument("-g", "--grid", dest="grid_spacing", type=int, help="horizontal grid resolution", default=1000)

options = parser.parse_args()
args = options.FILE
grid_spacing = options.grid_spacing

if len(args) == 0:
    nc_outfile = "melange_back_pressure_max.nc"
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

time = np.arange(0, 365) + 0.5


var = "time"
var_out = nc.createVariable(var, "d", dimensions=("time"))
var_out.axis = "T"
var_out.units = "day"
var_out.calendar = "365_day"
var_out.long_name = "time"
var_out.bounds = "time_bounds"
var_out[:] = time

var = "time_bounds"
var_out = nc.createVariable(var, "d", dimensions=("time", "nb"))
var_out.bounds = "time_bounds"
var_out[:, 0] = time - 0.5
var_out[:, 1] = time + 0.5


var = "frac_MBP"
var_out = nc.createVariable(var, "f", dimensions=("time"))
var_out.units = "1"

MBP_max = 1

winter_a = 300
winter_e = 90
spring_e = 105

MBP = np.zeros(len(time))
for k, t in enumerate(time):
    if t < 150:
        MBP[k] = MBP_max / np.sqrt(150) * np.sqrt(t)
    elif (t > 150) and (t < 200):
        MBP[k] = MBP_max - MBP_max / np.sqrt(200) * np.sqrt(t)
    else:
        MBP[k] = 0

var_out[:] = np.roll(MBP, -90)

nc.close()
