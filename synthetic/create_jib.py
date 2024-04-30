#!/usr/bin/env python
# Copyright (C) 2016-19 Andy Aschwanden

import numpy as np
from scipy.interpolate import griddata
from netCDF4 import Dataset as NC
from argparse import ArgumentParser


# set up the option parser
parser = ArgumentParser()
parser.description = "Generating synthetic outlet glacier / Jakobshavn Isbrae-like (adapted from Tinu Luethi)."
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
    nc_outfile = "jak" + str(grid_spacing) + "m.nc"
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

nmm_width = 10e3

H_m = 1200
x_t = 15e3
x_r = 7e3
beta_t = 2e-4
beta_r = -8e-5
mu = 5e-4

# Bed
Z_c = np.ones_like(X)

# Ellipsoid center
xe = x0
ye = 0
ze = 250

# Ellipsoid parameters
a = 500e3
b = 50e3
c = 250

Ze = -c * np.sqrt(
    1 - ((np.array(X, dtype=np.complex) - xe) / a) ** 2 - ((np.array(Y, dtype=np.complex) - ye) / b) ** 2
)
Ze = np.real(Ze) + ze

# Z_c = Ze

# Z_c[X < x_r] = 1 + beta_r * (x_r - X[X < x_r])
# Z_c[X >= x_t] = np.exp(-(beta_t * (x_t - X[X >= x_t])) ** 2)
# Z_b = -Z_c * np.exp(-(mu * Y) ** 2) * H_m

# Ellipsoid center
xc = 10e3
yc = 0
zc = 250

# Ellipsoid parameters
ac = 25e3
bc = 4.0e3
cc = 750


Zc = -cc * np.sqrt(
    1 - ((np.array(X, dtype=np.complex) - xc) / ac) ** 2 - ((np.array(Y, dtype=np.complex) - yc) / bc) ** 2
)
Zc = np.real(Zc) + zc

Z_b = Zc + Ze
# alpha = 0.10  # angle of plane
# Xp = np.cos(np.deg2rad(alpha)) * X - np.sin(np.deg2rad(alpha)) * Z_b
# Yp = Y
# Zp = np.sin(np.deg2rad(alpha)) * X + np.cos(np.deg2rad(alpha)) * Z_b
# # The only problem is now that Xp, Yp is no longer a regular grid, so you need to interpolate back onto the original grid:
# points = np.vstack(((np.ndarray.flatten(Xp), np.ndarray.flatten(Yp)))).T
# values = np.vstack((np.ndarray.flatten(Zp)))
# xi = np.vstack((np.ndarray.flatten(X), np.ndarray.flatten(Y))).T
# Zpi = griddata(points, values, xi, method="linear")

# Z_b = Zpi.reshape(N, M)
Z_b[X < 0] = -800

radius = 50e3
xcl, ycl = 25e3, y0 - 10e3
xcu, ycu = 25e3, y1 + 10e3
a = 125e3
b = 50e3

CL = (X - xcl) ** 2 / a ** 2 + (Y - ycl) ** 2 / b ** 2 < 1 ** 2
CU = (X - xcu) ** 2 / a ** 2 + (Y - ycu) ** 2 / b ** 2 < 1 ** 2

wall_elevation = 3000.0
if has_sidewalls:
    Z_b[np.logical_or(CL, CU)] = wall_elevation
    Z_b[np.logical_and((X < xcl), (Y < ycl + b))] = wall_elevation
    Z_b[np.logical_and((X < xcu), (Y > ycu - b))] = wall_elevation
    Z_b[X > x1 - nmm_width - m_buffer] = wall_elevation
    Z_b[Y < y0 + nmm_width + m_buffer] = wall_elevation
    Z_b[Y > y1 - nmm_width - m_buffer] = wall_elevation

# Surface
Z_s = np.zeros_like(X) + 0
Z_s[X >= 0] = 50 + 5.0 * np.sqrt(X[X >= 0])

thickness = Z_s - Z_b
thickness[thickness < 0] = 0
thickness[X < 0] = 0

# Create netCDF file
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
var_out = nc.createVariable(var, "f", dimensions=("y", "x"))
var_out.units = "meters"
var_out.standard_name = "bedrock_altitude"
var_out[:] = np.fliplr(Z_b[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx])

var = "usurf"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"))
var_out.units = "meters"
var_out.standard_name = "surface_altitude"
if has_sidewalls:
    Z_s[np.logical_or(CL, CU)] = wall_elevation
    Z_s[np.logical_and((X < xcl), (Y < ycl + radius))] = wall_elevation
    Z_s[np.logical_and((X < xcu), (Y > ycu - radius))] = wall_elevation

var_out[:] = np.fliplr(Z_s[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx])

var = "thk"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"))
var_out.units = "meters"
var_out.standard_name = "land_ice_thickness"
if has_sidewalls:
    thickness[np.logical_or(CL, CU)] = 0
    thickness[np.logical_and((X < xcl), (Y < ycl + radius))] = 0
    thickness[np.logical_and((X < xcu), (Y > ycu - radius))] = 0
    thickness[X < x0 + nmm_width + m_buffer] = 0
    thickness[X > x1 - nmm_width - m_buffer] = 0
    thickness[Y < y0 + nmm_width + m_buffer] = 0
    thickness[Y > y1 - nmm_width - m_buffer] = 0
var_out[:] = np.fliplr(thickness[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx])

var = "ftt_mask"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"))
var_out.units = ""
var_out.flag_meanings = "normal special_treatment"
var_out.long_name = "mask: zeros (modeling domain) and ones (no-model buffer near grid edges)"
var_out.flag_values = 0.0, 1.0
var_out.pism_intent = "model_state"
ftt_mask = np.zeros_like(Z_b)
ftt_mask[X < x0 + nmm_width + m_buffer] = 1
ftt_mask[X > x1 - nmm_width - m_buffer] = 1
ftt_mask[Y < y0 + nmm_width + m_buffer] = 1
ftt_mask[Y > y1 - nmm_width - m_buffer] = 1
if has_sidewalls:
    ftt_mask[np.logical_or(CL, CU)] = 1
    ftt_mask[np.logical_and((X < xcl), (Y < ycl + radius))] = 1
    ftt_mask[np.logical_and((X < xcu), (Y > ycu - radius))] = 1

mask = np.fliplr(ftt_mask.reshape(N, M)[m_buffer_idy:-m_buffer_idy, m_buffer_idx:-m_buffer_idx])
var_out[:] = mask

var = "no_model_mask"
var_out = nc.createVariable(var, "f", dimensions=("y", "x"))
var_out.units = ""
var_out.flag_meanings = "normal special_treatment"
var_out.long_name = "mask: zeros (modeling domain) and ones (no-model buffer near grid edges)"
var_out.flag_values = 0.0, 1.0
var_out.pism_intent = "model_state"

var_out[:] = mask


nc.close()
