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
    nc_outfile = "jak_runoff_" + str(grid_spacing) + "m.nc"
elif len(args) == 1:
    nc_outfile = args[0]
else:
    print("wrong number arguments, 0 or 1 arguments accepted")
    parser.print_help()
    import sys

    sys.exit(0)


dx = dy = grid_spacing  # m

# Domain extend
x0, x1 = -20.0e3, 150.0e3
y0, y1 = -25.0e3, 25.0e3

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
nc.createDimension("time", size=365)
nc.createDimension("nb", size=2)

t = np.arange(0, 365) + 0.5

runoff = np.zeros_like(X)

# runoff = np.sin(t / np.pi)
runoff = 1000 - 1000.0 / 100e3 * X
runoff[X < 0] = 0
runoff[X > 100e3] = 0

var = "time"
var_out = nc.createVariable(var, "d", dimensions=("time"))
var_out.axis = "T"
var_out.units = "days since 1-1-1"
var_out.calendar = "365_day"
var_out.long_name = "time"
var_out.bounds = "time_bounds"
var_out[:] = t

var = "time_bounds"
var_out = nc.createVariable(var, "d", dimensions=("time", "nb"))
var_out.bounds = "time_bounds"
var_out[:, 0] = t - 0.5
var_out[:, 1] = t + 0.5

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

for t in range(0, 365):
    var_out[t, :] = runoff[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx] * np.sin(
        (t + 0.5) * 1 * np.pi / 365
    )


nc.close()
