#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

from cftime import utime
from dateutil import rrule
from dateutil.parser import parse
from datetime import datetime

from netCDF4 import Dataset as NC
import pymc3 as pm
import numpy as np
import pandas as pd
import pylab as plt
from pyproj import Proj

from pymc3.gp.util import plot_gp_dist


def melting_point_temperature(depth, salinity):
    a = [-0.0575, 0.0901, -7.61e-4]
    return a[0] * salinity + a[1] + a[2] * depth


def create_nc(nc_outfile, theta_ocean, grid_spacing, time_dict):
    """
    Generate netCDF file
    """

    time = time_dict["time"]
    time_units = time_dict["units"]
    time_calendar = time_dict["calendar"]
    time_bnds = time_dict["time_bnds"]

    nt = len(theta_ocean)
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
    M = int((e1 - e0) / de) + 1
    N = int((n1 - n0) / dn) + 1

    easting = np.linspace(e0, e1, M)
    northing = np.linspace(n0, n1, N)
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
    gc_easting = np.zeros((M, grid_corners))
    # array holding y-component of grid corners
    gc_northing = np.zeros((N, grid_corners))
    # array holding the offsets from the cell centers
    # in x-direction (counter-clockwise)
    de_vec = np.array([-de / 2, de / 2, de / 2, -de / 2])
    # array holding the offsets from the cell centers
    # in y-direction (counter-clockwise)
    dn_vec = np.array([-dn / 2, -dn / 2, dn / 2, dn / 2])
    # array holding lat-component of grid corners
    gc_lat = np.zeros((N, M, grid_corners))
    # array holding lon-component of grid corners
    gc_lon = np.zeros((N, M, grid_corners))

    for corner in range(0, grid_corners):
        ## grid_corners in x-direction
        gc_easting[:, corner] = easting + de_vec[corner]
        # grid corners in y-direction
        gc_northing[:, corner] = northing + dn_vec[corner]
        # meshgrid of grid corners in x-y space
        gc_ee, gc_nn = np.meshgrid(gc_easting[:, corner], gc_northing[:, corner])
        # project grid corners from x-y to lat-lon space
        gc_lon[:, :, corner], gc_lat[:, :, corner] = proj(gc_ee, gc_nn, inverse=True)

    nc = NC(nc_outfile, "w", format="NETCDF4", compression_level=2)

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
    time_bnds_var = nc.createVariable(bnds_var_name, "d", dimensions=(time_dim, bnds_dim))
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

    var = "theta_ocean"
    var_out = nc.createVariable(var, "f", dimensions=("time", "y", "x"), fill_value=-2e9)
    var_out.units = "Celsius"
    var_out.long_name = "Just A Dummy"
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = np.repeat(theta_ocean, M * N).reshape(nt, N, M)

    mapping = nc.createVariable("mapping", "c")
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.0
    mapping.false_northing = 0.0
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = 90.0
    mapping.standard_parallel = 70.0
    mapping.straight_vertical_longitude_from_pole = -45.0

    # writing global attributes
    nc.Conventions = "CF 1.5"
    nc.close()


if __name__ == "__main__":

    depth = 250
    salinity = 34
    grid_spacing = 4500

    calendar = "standard"
    units = "days since 1980-1-1"
    cdftime_days = utime(units, calendar)

    start_date = datetime(1980, 1, 1)
    end_date = datetime(2020, 1, 1)
    end_date_yearly = datetime(2020, 1, 2)

    # create list with dates from start_date until end_date with
    # periodicity prule.
    bnds_datelist = list(rrule.rrule(rrule.MONTHLY, dtstart=start_date, until=end_date_yearly))
    bnds_datelist_yearly = list(rrule.rrule(rrule.YEARLY, dtstart=start_date, until=end_date_yearly))

    # calculate the days since refdate, including refdate, with time being the
    bnds_interval_since_refdate = cdftime_days.date2num(bnds_datelist)
    bnds_interval_since_refdate_yearly = cdftime_days.date2num(bnds_datelist_yearly)
    time_interval_since_refdate = bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2

    time_dict = {
        "calendar": calendar,
        "units": units,
        "time": time_interval_since_refdate,
        "time_bnds": bnds_interval_since_refdate,
    }

    dpy = np.diff(bnds_interval_since_refdate_yearly)
    dpy = np.repeat(dpy, 12, axis=0)

    decimal_time = bnds_interval_since_refdate[0:-1] / dpy + 1980

    kh_df = pd.read_csv("disko-bay_khazendar.csv", names=["time", "temperature"])
    mo_df = pd.read_csv("disko-bay_motyka.csv", names=["time", "temperature"])
    all_df = pd.concat([kh_df, mo_df])
    all_df = all_df.sort_values(by="time")

    X_kh = kh_df.time.values.reshape(-1, 1)
    X_mo = mo_df.time.values.reshape(-1, 1)
    y_kh = kh_df.temperature.values
    y_mo = mo_df.temperature.values

    X = all_df.time.values.reshape(-1, 1)
    y = all_df.temperature.values

    X_new = decimal_time[:, None]

    f = "f_pred"

    covs = {
        "mat32": pm.gp.cov.Matern32,
        "mat52": pm.gp.cov.Matern52,
        "exp": pm.gp.cov.Exponential,
        "exp-quad": pm.gp.cov.ExpQuad,
    }

    for kernel, cov in covs.items():
        print(kernel)
        # Khazendar and Motyka merged
        with pm.Model() as gp:
            ρ = pm.Normal("ρ", 1)
            η = pm.Normal("η", 5)
            K = η * cov(1, ρ)
            σ = pm.Normal("σ", 0.2)
            gp = pm.gp.Marginal(cov_func=K)
            y_ = gp.marginal_likelihood("y", X=X, y=y, noise=σ)
            mp = pm.find_MAP()
            f_pred = gp.conditional(f, X_new)
            pred_samples = pm.sample_posterior_predictive([mp], vars=[f_pred], samples=10)

        for s, temperate in enumerate(pred_samples[f]):
            theta_ocean = temperate - melting_point_temperature(depth, salinity)
            ofile = f"disko_bay_theta_ocean_{kernel}_{s}_1980_2019.nc"
            create_nc(ofile, theta_ocean, grid_spacing, time_dict)

        fig = plt.figure()
        ax = fig.gca()

        # plot the samples from the gp posterior with samples and shading
        plot_gp_dist(ax, pred_samples[f], X_new, palette="Blues")

        # plot the data and the true latent function
        ax.plot(X_kh, y_kh, "o", color="0.4", ms=5, label="Observed data Khazendar")
        ax.plot(X_mo, y_mo, "o", color="0.6", ms=5, label="Observed data Motyka")

        ax.set_xlabel("Time")
        ax.set_ylabel("Temperature (Celsius)")
        ax.set_xlim(1980, 2020)
        ax.set_ylim(0)
        plt.legend()
        fig.savefig(f"disko-bay-temps_{kernel}.pdf")
