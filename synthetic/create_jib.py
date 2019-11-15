#!/usr/bin/env python
# Copyright (C) 2016, 2017 Andy Aschwanden

import numpy as np
from scipy.interpolate import griddata
from netCDF4 import Dataset as NC
from argparse import ArgumentParser


# set up the option parser
parser = ArgumentParser()
parser.description = "Generating synthetic outlet glacier."
parser.add_argument("FILE", nargs="*")
parser.add_argument("-g", "--grid", dest="grid_spacing", type=int, help="horizontal grid resolution", default=1000)
parser.add_argument(
    "-s", "--side_walls", dest="has_sidewalls", action="store_true", help="horizontal grid resolution", default=False
)

options = parser.parse_args()
args = options.FILE
grid_spacing = options.grid_spacing
has_sidewalls = options.has_sidewalls

if len(args) == 0:
    nc_outfile = "og" + str(grid_spacing) + "m.nc"
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
y0, y1 = -50.0e3, 50.0e3

# shift to cell centers
x0 += dx / 2
x1 -= dx / 2
y0 += dy / 2
y1 -= dy / 2

# Number of grid points
M = int((x1 - x0) / dx) + 1
N = int((y1 - y0) / dy) + 1

x = np.linspace(x0, x1, M)
y = np.linspace(y0, y1, N)
X, Y = np.meshgrid(x, y)

H_m = 1600
x_t = 15
x_r = 7
beta_t = 2e-4
beta_r = -8e-5
mu = 7e-4

z_c = np.ones_like(X)

z_c[X < x_r] = (1 + beta_r) * (x_r - X[X < x_r])
z_c[X >= x_t] = np.exp(-(beta_t * (x_t - X[X >= x_t])) ** 2)
z_b = -z_c * np.exp(-(mu * Y) ** 2) * H_m

z_s = np.zeros_like(X) + 50
z_s[X >= 0] = 50 + 4.5 * np.sqrt(X[X >= 0])

nc = NC(nc_outfile, "w", format="NETCDF4")

nc.createDimension("x", size=M)
nc.createDimension("y", size=N)

var = "x"
var_out = nc.createVariable(var, "d", dimensions=("x"))
var_out.axis = "X"
var_out.long_name = "X-coordinate in Cartesian system"
var_out.standard_name = "projection_x_coordinate"
var_out.units = "meters"
var_out[:] = x

var = "y"
var_out = nc.createVariable(var, "d", dimensions=("y"))
var_out.axis = "Y"
var_out.long_name = "Y-coordinate in Cartesian system"
var_out.standard_name = "projection_y_coordinate"
var_out.units = "meters"
var_out[:] = y

var = "topg"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"))
var_out.units = "meters"
var_out.standard_name = "bedrock_altitude"
var_out[:] = z_b

var = "usurf"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"), fill_value=0)
var_out.units = "meters"
var_out.standard_name = "surface_altitude"
var_out[:] = z_s
var = "thk"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"), fill_value=0)
var_out.units = "meters"
var_out.standard_name = "land_ice_thickness"
var_out[:] = z_s - z_b


nc.close()
