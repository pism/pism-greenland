#!/bin/env python3
# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
from netCDF4 import Dataset as NC
import numpy as np

# set up the option parser
parser = ArgumentParser()
parser.add_argument("FILE", nargs=1)

options = parser.parse_args()
m_file = options.FILE[0]

nc = NC(m_file, "a")
x = nc.variables["x"][:]
y = nc.variables["y"][:]
nx = len(x)
ny = len(y)

X, Y = np.meshgrid(x, y)

A = 500.0
sigma = 1000.0
x0, y0 = -197000.0, -2275300.0

Z = -520 + A * np.exp(-((X - x0) ** 2) / (2 * sigma ** 2) - ((Y - y0) ** 2) / (2 * sigma ** 2))
mask = np.zeros((ny, nx))
mask[(X - x0) ** 2 + (Y - y0) ** 2 < A ** 2] = 1
bed_elevation = nc.variables["bed"][:]
bed_elevation = np.where(mask, Z, bed_elevation)
nc.variables["bed"][:] = bed_elevation
nc.close()
