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
import time
from itertools import product
import os
import statsmodels.api as sm


def toDecimalYear(date):
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())

    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year + 1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return date.year + fraction


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

    # depths to average over
    depth_min = 225
    depth_max = 275
    # depth for freezing point calculation
    depth = 250
    salinity = 34
    grid_spacing = 4500

    calendar = "standard"
    units = "days since 1980-1-1"
    cdftime_days = utime(units, calendar)

    start_date = datetime(1980, 1, 1)
    end_date = datetime(2021, 1, 1)
    end_date_yearly = datetime(2021, 1, 2)

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

    mo_df = pd.read_csv("disko_bay_motyka.csv")

    omg_df = pd.read_csv("disko_bay_omg_axctd.csv", na_values=-99.0).dropna()
    omg_df = omg_df[(omg_df["Depth"] <= depth_max) & (omg_df["Depth"] >= depth_min)]
    omg_time = pd.to_datetime(omg_df.Date, format="%m/%d/%Y %H:%M:%S")
    omg_df["Date"] = omg_time
    omg_df = omg_df.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_df = omg_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    omg_time = [toDecimalYear(d) for d in omg_df.index]
    omg_df["Date"] = omg_time

    ices_df = pd.read_csv("disko_bay_ices.csv")
    ices_df = ices_df[(ices_df.Depth >= depth_min) & (ices_df.Depth <= depth_max)].reset_index(drop=True)
    ices_time = pd.to_datetime(ices_df.Date, format="%Y/%m/%d %H:%M:%S")
    ices_df["Date"] = ices_time
    ices_df = ices_df.set_index("Date")
    ices_df = ices_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    ices_time = [toDecimalYear(d) for d in ices_df.index]
    ices_df["Date"] = ices_time
    ices_df.to_csv("disko_bay_ices_depth_averaged.csv")

    holl_df = pd.read_csv("disko_bay_xctd_holland.csv", parse_dates=[0])
    holl_df = holl_df[(holl_df.Depth >= depth_min) & (holl_df.Depth <= depth_max)].reset_index(drop=True)
    holl_df = holl_df.set_index("Date")
    holl_df = holl_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    holl_time = [toDecimalYear(d) for d in holl_df.index]
    holl_df["Date"] = holl_time
    holl_df.to_csv("disko_bay_ices_depth_averaged.csv")

    all_df = pd.concat([mo_df, ices_df, omg_df, holl_df])
    all_df = all_df.sort_values(by="Date")

    X_mo = mo_df.Date.values.reshape(-1, 1)
    X_ices = ices_df.Date.values.reshape(-1, 1)
    X_omg = omg_df.Date.values.reshape(-1, 1)
    X_holl = holl_df.Date.values.reshape(-1, 1)

    y_mo = mo_df.Temperature.values
    y_ices = ices_df.Temperature.values
    y_omg = omg_df.Temperature.values
    y_holl = holl_df.Temperature.values

    trend_start = [1980, 1998]
    trend_end = [1998, 2016]

    df = all_df
    x_var = "Date"
    y_var = "Temperature"

    fig = plt.figure()
    ax = fig.gca()

    ax.plot(X_holl, y_holl, "o", color="#f4a582", ms=4, label="Observed (Holland)")
    ax.plot(X_mo, y_mo, "o", color="#b2182b", ms=4, label="Observed (Motyka)")
    ax.plot(X_ices, y_ices, "o", color="#92c5de", ms=4, label="Observed (ICES)")
    ax.plot(X_omg, y_omg, "o", color="#2166ac", ms=4, label="Observed (OMG)")

    # for a, e in zip(trend_start, trend_end):
    #     m_df = df[(df[x_var] >= a) & (df[x_var] <= e)]
    #     x = m_df[x_var]
    #     y = m_df[y_var]
    #     X = sm.add_constant(x)
    #     ols = sm.OLS(y, X).fit()
    #     p = ols.params
    #     bias = p[0]
    #     trend = p[1]
    #     trend_error = ols.bse[-1]
    #     mean = np.mean(y)
    #     ax.plot([a, e], np.array([a, e]) * trend + bias, linewidth=2, color="k")
    #     ax.plot([a, e], np.array([mean, mean]), linewidth=2, color="0.5")
        
    ax.set_xlabel("Time")
    ax.set_ylabel("Temperature (Celsius)")
    ax.set_xlim(1980, 2021)
    ax.set_ylim(0, 5)
    plt.legend()
    fig.savefig("disko-bay-temps-obs.pdf")
    plt.cla()
    plt.close(fig)

    X = all_df.Date.values.reshape(-1, 1)
    y = all_df.Temperature.values

    # Normalize
    X_mean = X.mean(axis=0)
    y_mean = y.mean(axis=0)
    # X -= X_mean
    # y -= y_mean

    X_new = decimal_time[:, None]

    f = "f_pred"

    # covs = {
    #     "mat32": pm.gp.cov.Matern32,
    #     "mat52": pm.gp.cov.Matern52,
    #     "exp": pm.gp.cov.Exponential,
    #     "exp-quad": pm.gp.cov.ExpQuad,
    # }

    covs = {
        "mat32": pm.gp.cov.Matern32,
    }


    odir = "test"
    if not os.path.isdir(odir):
        os.makedirs(odir)
    mus = [0.2, 1.00]
    sigmas = [0.2, 1.00]
    combinations = list(product(mus, sigmas, mus, sigmas))
#    combinations = [(2, 1, 1, 1)]
    for kernel, cov in covs.items():
        print(f"Covariance function: {kernel}")
        for rho_mu, rho_sigma, eta_mu, eta_sigma in combinations:
            print(f"rho_mu={rho_mu:.2f}, rho_sigma={rho_sigma:.2f} eta_mu={eta_mu:.2f}, eta_sigma={eta_sigma:.2f}")
            with pm.Model() as gp:
                ℓ = pm.Normal("ℓ", mu=rho_mu, sigma=rho_sigma)
                η = pm.Normal("η", mu=eta_mu, sigma=eta_sigma)
                mean_func = pm.gp.mean.Zero()
                cov_func = η * cov(1, ℓ)
                noise_mu = 0.05
                noise_sigma = 0.1
                σ = pm.Normal("σ", sigma=noise_sigma, mu=noise_mu)
                gp = pm.gp.Marginal(mean_func=mean_func, cov_func=cov_func)
                y_ = gp.marginal_likelihood("y", X=X, y=y, noise=σ)
                mp = pm.find_MAP()
                mp_df = pd.DataFrame(
                    {
                        "Parameter": ["ℓ", "η", "σ"],
                        "Value at MAP": [float(mp["ℓ"]), float(mp["η"]), float(mp["σ"])],
                    }
                )
                print(mp_df)
                f_pred = gp.conditional("f_pred", X_new, pred_noise=False)
                f_samples = pm.sample_posterior_predictive([mp], var_names=["f_pred"], samples=10)
                y_pred = gp.conditional("y_pred", X_new, pred_noise=True)
                y_samples = pm.sample_posterior_predictive([mp], var_names=["y_pred"], samples=10)

            # for s, temperate in enumerate(pred_samples[f]):
            #     theta_ocean = temperate - melting_point_temperature(depth, salinity)
            #     ofile = f"disko_bay_theta_ocean_{kernel}_{s}_1980_2019.nc"
            #     create_nc(ofile, theta_ocean, grid_spacing, time_dict)

            fig = plt.figure()
            ax = fig.gca()

            # plot the samples from the gp posterior with samples and shading
            plot_gp_dist(ax, f_samples["f_pred"], X_new, palette="Greys")

            # plot the data and the true latent function
            ax.plot(X_holl, y_holl, "o", color="#f4a582", ms=4, label="Observed (Holland)")
            ax.plot(X_mo, y_mo, "o", color="#b2182b", ms=4, label="Observed (Motyka)")
            ax.plot(X_ices, y_ices, "o", color="#92c5de", ms=4, label="Observed (ICES)")
            ax.plot(X_omg, y_omg, "o", color="#2166ac", ms=4, label="Observed (OMG)")

            ax.set_xlabel("Time")
            ax.set_ylabel("Temperature (Celsius)")
            ax.set_xlim(1980, 2021)
            ax.set_ylim(0, 5)
            plt.legend()
            fig.savefig(
                f"{odir}/disko-bay-temps_{kernel}_rho_mu_{rho_mu:.2f}_rho_sigma_{rho_sigma:.2f}_eta_mu_{eta_mu:.2f}_eta_sigma_{eta_sigma:.2f}_noise_mu_{noise_mu:.2f}_noise_sigma_{noise_sigma:.2f}_no_noise.pdf"
            )
            plt.cla()
            plt.close(fig)

            fig = plt.figure()
            ax = fig.gca()

            # plot the samples from the gp posterior with samples and shading
            plot_gp_dist(ax, y_samples["y_pred"], X_new, palette="Greys")

            # plot the data and the true latent function
            ax.plot(X_holl, y_holl, "o", color="#f4a582", ms=4, label="Observed (Holland)")
            ax.plot(X_mo, y_mo, "o", color="#b2182b", ms=4, label="Observed (Motyka)")
            ax.plot(X_ices, y_ices, "o", color="#92c5de", ms=4, label="Observed (ICES)")
            ax.plot(X_omg, y_omg, "o", color="#2166ac", ms=4, label="Observed (OMG)")

            ax.set_xlabel("Time")
            ax.set_ylabel("Temperature (Celsius)")
            ax.set_xlim(1980, 2021)
            ax.set_ylim(0, 5)
            plt.legend()
            fig.savefig(
                f"{odir}/disko-bay-temps_{kernel}_rho_mu_{rho_mu:.2f}_rho_sigma_{rho_sigma:.2f}_eta_mu_{eta_mu:.2f}_eta_sigma_{eta_sigma:.2f}_noise_mu_{noise_mu:.2f}_noise_sigma_{noise_sigma:.2f}_w_noise.pdf"
            )
            plt.cla()
            plt.close(fig)
