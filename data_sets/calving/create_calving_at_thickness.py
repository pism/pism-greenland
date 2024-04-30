#!/usr/bin/env python3
# Copyright (C) 2020-21 Andy Aschwanden
#
# Generate time-series of calving in Jakobshavn Fjord

from argparse import ArgumentParser

import cftime
from dateutil import rrule
from datetime import datetime
from netCDF4 import Dataset as NC
import numpy as np
import pandas as pd
from pyproj import Proj




if __name__ == "__main__":
    __spec__ = None

    # set up the option parser
    parser = ArgumentParser()
    parser.add_argument("FILE", nargs=1)
    parser.add_argument(
        "--tct_0", dest="tct_0", type=float, help="southern thickness calving threshold, in m", default=400
    )
    parser.add_argument("--tct_1", dest="tct_1", type=float, help="northern thickness calving threshold, in m", default=50)
    parser.add_argument(
        "--lat_0", dest="lat_0", type=float, help="latitude to apply southern thickness calving threshold", default=74
    )
    parser.add_argument(
        "--lat_1", dest="lat_1", type=float, help="latitude to apply northern thickness calving threshold", default=76
    )

    options = parser.parse_args()
    args = options.FILE
    ofile = args[0]
    lat_0 = options.lat_0
    lat_1 = options.lat_1
    tct_0 = options.tct_0
    tct_1 = options.tct_1
    
    a_tct = (tct_1 - tct_0) / (lat_1 - lat_0)
    b_tct = tct_0 - a_tct * lat_0

    grid_spacing = 1200

    # The temporal averaging window
    freq = "1D"

    start_date = datetime(1980, 1, 1)
    end_date = datetime(2021, 1, 1)
    end_date_yearly = datetime(2021, 1, 2)

    calendar = "standard"
    units = "days since 1980-1-1"

    # create list with dates from start_date until end_date with
    # periodicity prule for netCDF file
    # and use data_range for pytorch X_new.
    sampling_interval = "yearly"
    dates = pd.date_range(start=start_date, end=end_date, freq="1AS")

    rd = {
        "daily": rrule.DAILY,
        "weekly": rrule.WEEKLY,
        "monthly": rrule.MONTHLY,
        "yearly": rrule.YEARLY,
    }

    bnds_datelist = list(
        rrule.rrule(rd[sampling_interval], dtstart=start_date, until=end_date_yearly)
    )
    # calculate the days since refdate, including refdate, with time being the
    bnds_interval_since_refdate = cftime.date2num(
        bnds_datelist, units, calendar=calendar
    )
    time_interval_since_refdate = (
        bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2
    )

    time_dict = {
        "calendar": calendar,
        "units": units,
        "time": time_interval_since_refdate,
        "time_bnds": bnds_interval_since_refdate,
    }

    """
    Generate netCDF file
    """

    time_units = time_dict["units"]
    time_calendar = time_dict["calendar"]
    time = time_dict["time"]
    time_bnds = time_dict["time_bnds"]

    nt = len(dates)
    xdim = "x"
    ydim = "y"

    # define output grid, these are the extents of Mathieu's domain (cell
    # corners)
    e0 = -638000
    n0 = -3349600
    e1 = 864700
    n1 = -657600

    # Add a buffer on each side such that we get nice grids up to a grid spacing
    # of 36 km.

    buffer_e = 148650
    buffer_n = 130000
    e0 -= buffer_e + 468000
    n0 -= buffer_n
    e1 += buffer_e
    n1 += buffer_n

    # Shift to cell centers
    e0 += grid_spacing / 2
    n0 += grid_spacing / 2
    e1 -= grid_spacing / 2
    n1 -= grid_spacing / 2

    de = dn = grid_spacing  # m
    m = int((e1 - e0) / de) + 1
    n = int((n1 - n0) / dn) + 1

    easting = np.linspace(e0, e1, m)
    northing = np.linspace(n0, n1, n)
    ee, nn = np.meshgrid(easting, northing)

    # Set up EPSG 3413 (NSIDC north polar stereo) projection
    projection = "epsg:3413"
    proj = Proj(projection)

    lon, lat = proj(ee, nn, inverse=True)

    # number of grid corners
    grid_corners = 4
    # grid corner dimension name
    grid_corner_dim_name = "nv4"

    # array holding x-component of grid corners
    gc_easting = np.zeros((m, grid_corners))
    # array holding y-component of grid corners
    gc_northing = np.zeros((n, grid_corners))
    # array holding the offsets from the cell centers
    # in x-direction (counter-clockwise)
    de_vec = np.array([-de / 2, de / 2, de / 2, -de / 2])
    # array holding the offsets from the cell centers
    # in y-direction (counter-clockwise)
    dn_vec = np.array([-dn / 2, -dn / 2, dn / 2, dn / 2])
    # array holding lat-component of grid corners
    gc_lat = np.zeros((n, m, grid_corners))
    # array holding lon-component of grid corners
    gc_lon = np.zeros((n, m, grid_corners))

    for corner in range(0, grid_corners):
        # grid_corners in x-direction
        gc_easting[:, corner] = easting + de_vec[corner]
        # grid corners in y-direction
        gc_northing[:, corner] = northing + dn_vec[corner]
        # meshgrid of grid corners in x-y space
        gc_ee, gc_nn = np.meshgrid(gc_easting[:, corner], gc_northing[:, corner])
        # project grid corners from x-y to lat-lon space
        gc_lon[:, :, corner], gc_lat[:, :, corner] = proj(gc_ee, gc_nn, inverse=True)

    nc = NC(ofile, "w", format="NETCDF4", compression_level=2)

    nc.createDimension(xdim, size=easting.shape[0])
    nc.createDimension(ydim, size=northing.shape[0])

    time_dim = "time"
    if time_dim not in list(nc.dimensions.keys()):
        nc.createDimension(time_dim)

    # create a new dimension for bounds only if it does not yet exist
    bnds_dim = "nb2"
    if bnds_dim not in list(nc.dimensions.keys()):
        nc.createDimension(bnds_dim, 2)

    # variable names consistent with PISM
    time_var_name = "time"
    bnds_var_name = "time_bnds"

    # create time variable
    time_var = nc.createVariable(time_var_name, "d", dimensions=(time_dim))
    time_var[:] = time
    time_var.bounds = bnds_var_name
    time_var.units = time_units
    time_var.calendar = time_calendar
    time_var.standard_name = time_var_name
    time_var.axis = "T"

    # create time bounds variable
    time_bnds_var = nc.createVariable(
        bnds_var_name, "d", dimensions=(time_dim, bnds_dim)
    )
    time_bnds_var[:, 0] = time_bnds[0:-1]
    time_bnds_var[:, 1] = time_bnds[1::]

    var = xdim
    var_out = nc.createVariable(var, "d", dimensions=(xdim))
    var_out.axis = xdim
    var_out.long_name = "X-coordinate in Cartesian system"
    var_out.standard_name = "projection_x_coordinate"
    var_out.units = "meters"
    var_out[:] = easting

    var = ydim
    var_out = nc.createVariable(var, "d", dimensions=(ydim))
    var_out.axis = ydim
    var_out.long_name = "Y-coordinate in Cartesian system"
    var_out.standard_name = "projection_y_coordinate"
    var_out.units = "meters"
    var_out[:] = northing

    var = "lon"
    var_out = nc.createVariable(var, "d", dimensions=(ydim, xdim))
    var_out.units = "degrees_east"
    var_out.valid_range = -180.0, 180.0
    var_out.standard_name = "longitude"
    var_out.bounds = "lon_bnds"
    var_out[:] = lon

    var = "lat"
    var_out = nc.createVariable(var, "d", dimensions=(ydim, xdim))
    var_out.units = "degrees_north"
    var_out.valid_range = -90.0, 90.0
    var_out.standard_name = "latitude"
    var_out.bounds = "lat_bnds"
    var_out[:] = lat

    nc.createDimension(grid_corner_dim_name, size=grid_corners)

    var = "lon_bnds"
    # Create variable 'lon_bnds'
    var_out = nc.createVariable(var, "f", dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lon_bnds'
    var_out.units = "degreesE"
    # Assign values to variable 'lon_nds'
    var_out[:] = gc_lon

    var = "lat_bnds"
    # Create variable 'lat_bnds'
    var_out = nc.createVariable(var, "f", dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lat_bnds'
    var_out.units = "degreesN"
    # Assign values to variable 'lat_bnds'
    var_out[:] = gc_lat

    tct = a_tct * lat + b_tct
    tct[lat < lat_0] = a_tct * lat_0 + b_tct
    tct[lat > lat_1] = a_tct * lat_1 + b_tct

    
    var = "thickness_calving_threshold"
    var_out = nc.createVariable(
        var, "f", dimensions=("time", "y", "x"), fill_value=-2e9, zlib=True, complevel=2
    )
    var_out.units = "m"
    var_out.long_name = "threshold used by the 'calving at threshold' calving method"
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    for k, date in enumerate(dates):
        print(date)
        var_out[k, :] = tct

    mapping = nc.createVariable("mapping", "c")
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.0
    mapping.false_northing = 0.0
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = 90.0
    mapping.standard_parallel = 70.0
    mapping.straight_vertical_longitude_from_pole = -45.0

    # writing global attributes
    nc.Conventions = "CF-1.7"
    nc.close()
