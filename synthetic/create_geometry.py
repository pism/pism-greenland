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
x0, x1 = 0.0, 360.0e3
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

# Ellipsoid center
xe = 150e3 + m_buffer
ye = 0
ze = -750

# Ellipsoid parameters
a = 100e3
b = 7.5e3
c = 750

Ze = -c * np.sqrt(
    1 - ((np.array(X, dtype=np.complex) - xe) / a) ** 2 - ((np.array(Y, dtype=np.complex) - ye) / b) ** 2
)
Ze = np.real(Ze) + ze

# That works because outside the area of the ellipsoid the sqrt is purely imaginary (hence the 'real' command).
# Rotating this around the y-axis is not really complicated, per se. Rotation around an angle alpha would give you the following:

alpha = 0.25  # angle of ellipsoid
beta = 0.40  # angle of surface plane
Xp = np.cos(np.deg2rad(alpha)) * X - np.sin(np.deg2rad(alpha)) * Ze
Yp = Y
Zp = np.sin(np.deg2rad(alpha)) * X + np.cos(np.deg2rad(alpha)) * Ze
# The only problem is now that Xp, Yp is no longer a regular grid, so you need to interpolate back onto the original grid:
points = np.vstack(((np.ndarray.flatten(Xp), np.ndarray.flatten(Yp)))).T
values = np.vstack((np.ndarray.flatten(Zp)))
xi = np.vstack((np.ndarray.flatten(X), np.ndarray.flatten(Y))).T
Zpi = griddata(points, values, xi, method="linear")

a_s = 20
x_s = 60e3
Zs = np.real(np.sqrt(a_s * (np.array(X - x_s, dtype=np.complex))))
Xsp = np.cos(np.deg2rad(alpha)) * (X - x_s) - np.sin(np.deg2rad(alpha)) * Zs
Ysp = Y
Zsp = np.sin(np.deg2rad(alpha)) * (X - x_s) + np.cos(np.deg2rad(alpha)) * Zs
# The only problem is now that Xp, Yp is no longer a regular grid, so you need to interpolate back onto the original grid:
points = np.vstack(((np.ndarray.flatten(Xsp), np.ndarray.flatten(Ysp)))).T
values = np.vstack((np.ndarray.flatten(Zsp)))
xi = np.vstack((np.ndarray.flatten(X - x_s), np.ndarray.flatten(Y))).T

Zspi = griddata(points, values, xi, method="linear") - 250


radius = 25e3
xcl, ycl = 50e3, y0
xcu, ycu = 50e3, y1

CL = (X - xcl) ** 2 + (Y - ycl) ** 2 < radius ** 2
CU = (X - xcu) ** 2 + (Y - ycu) ** 2 < radius ** 2

wall_elevation = 1000.0
if has_sidewalls:
    Zpi[np.logical_or(CL, CU)] = wall_elevation
    Zpi[np.logical_and((X < xcl), (Y < ycl + radius))] = wall_elevation
    Zpi[np.logical_and((X < xcu), (Y > ycu - radius))] = wall_elevation

mZ = (Zspi - Zpi).T
thk = np.reshape(Zspi - Zpi, (N, M))

thk[X < x_s] = 0.0
thk[thk < 0] = 0.0

nc = NC(nc_outfile, "w", format="NETCDF4")

nc.createDimension("x", size=MM)
nc.createDimension("y", size=NN)

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

var = "topg"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"), fill_value=0)
var_out.units = "meters"
var_out.standard_name = "bedrock_altitude"
var_out[:] = Zpi.reshape(N, M)[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx]

var = "usurf"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"), fill_value=0)
var_out.units = "meters"
var_out.standard_name = "surface_altitude"
var_out[:] = Zspi.reshape(N, M)[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx]

var = "thk"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"), fill_value=0)
var_out.units = "meters"
var_out.standard_name = "land_ice_thickness"
var_out[:] = thk.reshape(N, M)[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx]

var = "no_model_mask"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"))
var_out.units = ""
var_out.flag_meanings = "normal special_treatment"
var_out.long_name = "mask: zeros (modeling domain) and ones (no-model buffer near grid edges)"
var_out.flag_values = 0.0, 1.0
var_out.pism_intent = "model_state"
no_model_mask = np.zeros_like(thk)
no_model_mask[X >= 256e3] = 1
var_out[:] = no_model_mask.reshape(N, M)[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx]

nc.close()
