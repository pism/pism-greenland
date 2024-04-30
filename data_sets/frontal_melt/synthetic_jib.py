#!/usr/bin/env python
# Copyright (C) 2019 Andy Aschwanden

import numpy as np
from scipy.interpolate import griddata
from netCDF4 import Dataset as NC
from argparse import ArgumentParser


# set up the option parser
parser = ArgumentParser()
parser.description = "Generating synthetic outlet glacier."
parser.add_argument("FILE", nargs="*")
parser.add_argument("-g", "--grid", dest="grid_spacing", type=int, help="horizontal grid resolution", default=1000)

options = parser.parse_args()
args = options.FILE
grid_spacing = options.grid_spacing

if len(args) == 0:
    nc_outfile = "synth_jib_runoff_g" + str(grid_spacing) + "m.nc"
elif len(args) == 1:
    nc_outfile = args[0]
else:
    print("wrong number arguments, 0 or 1 arguments accepted")
    parser.print_help()
    import sys

    sys.exit(0)


dx = dy = grid_spacing  # m

# Domain extend
x0, x1 = -50.0e3, 250.0e3
y0, y1 = -50.0e3, 50.0e3

# shift to cell centers
x0 += dx / 2
x1 -= dx / 2
y0 += dy / 2
y1 -= dy / 2

m_buffer = 10e3
# add buffer
x0 -= m_buffer
y0 -= m_buffer
x1 += m_buffer
y1 += m_buffer

# Number of grid points
M = int((x1 - x0) / dx) + 1
N = int((y1 - y0) / dy) + 1

x = np.linspace(x0, x1, M)
y = np.linspace(y0, y1, N)
X, Y = np.meshgrid(x, y)

m_buffer_idx = int(m_buffer / dx)
m_buffer_idy = int(m_buffer / dy)

MM = int(M - m_buffer / dx * 2)
NN = int(N - m_buffer / dy * 2)

# Create netCDF file
nc = NC(nc_outfile, "w", format="NETCDF4")

nc.createDimension("x", size=MM)
nc.createDimension("y", size=NN)
nc.createDimension("time")
nc.createDimension("nb", size=2)

time = np.arange(0, 3650) + 0.5

m_runoff = 4.0 * 910.0  # kg m-2 year-1
x_r = 500e3  # m

runoff = m_runoff - m_runoff / x_r * X
runoff[X < 0] = 0
runoff[X > x_r] = 0
runoff = np.fliplr(runoff)

m_theta = 4.0
theta = np.fliplr(np.zeros_like(X) + m_theta)

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

var = "x"
var_out = nc.createVariable(var, "d", dimensions=("x"))
var_out.axis = "X"
var_out.long_name = "X-coordinate in Cartesian system"
var_out.standard_name = "projection_x_coordinate"
var_out.units = "meters"
var_out[:] = x[m_buffer_idx:-m_buffer_idx]

var = "y"
var_out = nc.createVariable(var, "d", dimensions=("y"))
var_out.axis = "Y"
var_out.long_name = "Y-coordinate in Cartesian system"
var_out.standard_name = "projection_y_coordinate"
var_out.units = "meters"
var_out[:] = y[m_buffer_idy:-m_buffer_idy]

var = "water_input_rate"
var_out = nc.createVariable(var, "f", dimensions=("time", "y", "x"))
var_out.units = "kg m-2 year-1"

for year in range(0, 10):
    for day in range(0, 365):
        var_out[year * 365 + day, :] = 0
        if (day > 120) and (day < 245):
            var_out[year * 365 + day, :] = runoff[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx]

var = "theta_ocean"
var_out = nc.createVariable(var, "f", dimensions=("time", "y", "x"))
var_out.units = "Celsius"

for t in time:
    var_out[t, :] = theta[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx]

var = "salinity_ocean"
var_out = nc.createVariable(var, "f", dimensions=("time", "y", "x"))
var_out.units = "g/kg"
var_out = 34.0

nc.close()
