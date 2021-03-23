#!/usr/bin/env python
# Copyright (C) 2020-21 Andy Aschwanden

from cftime import utime
from dateutil import rrule
from datetime import datetime
from netCDF4 import Dataset as NC
import gpytorch
import torch
import numpy as np
import pandas as pd
import pylab as plt
from pyproj import Proj


class MultitaskGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, num_tasks):
        super(MultitaskGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_modules = [gpytorch.means.ConstantMean() for i in range(num_tasks)]
        self.covar_module = gpytorch.kernels.RBFKernel()

        # Surprisingly the Gram matrix of a rank-1 outer product appears to be sufficient
        # for parameterizing the inter-task covariance matrix, as increasing the rank
        # does not improved the fit.
        self.task_covar_module = gpytorch.kernels.IndexKernel(num_tasks=num_tasks, rank=1)

    def forward(self, x, i):
        # This is a hack that allows for different means to be queried when there are different task indices
        mean_x = torch.cat([self.mean_modules[ii](xx) for ii, xx in zip(i, x)])

        # Get input-input covariance
        covar_x = self.covar_module(x)
        # Get task-task covariance
        covar_i = self.task_covar_module(i)
        # Multiply the two together to get the covariance we want
        covar = covar_x.mul(covar_i)

        return gpytorch.distributions.MultivariateNormal(mean_x, covar)


def set_size(w, h, ax=None):
    """ w, h: width, height in inches """

    if not ax:
        ax = plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w) / (r - l)
    figh = float(h) / (t - b)
    ax.figure.set_size_inches(figw, figh)


def to_decimal_year(date):
    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction


def melting_point_temperature(depth, salinity):
    a = [-0.0575, 0.0901, -7.61e-4]
    return a[0] * salinity + a[1] + a[2] * depth


def create_nc(nc_outfile, theta, salinity, grid_spacing, start_date, end_date, calendar, untis):
    """
    Generate netCDF file
    """

    cdftime_days = utime(units, calendar)

    # create list with dates from start_date until end_date with
    # periodicity prule.
    bnds_datelist = list(rrule.rrule(rrule.MONTHLY, dtstart=start_date, until=end_date_yearly))

    # calculate the days since refdate, including refdate, with time being the
    bnds_interval_since_refdate = cdftime_days.date2num(bnds_datelist)
    time_interval_since_refdate = bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2

    time_units = units
    time_calendar = calendar
    time = time_interval_since_refdate
    time_bnds = bnds_interval_since_refdate

    nt = len(theta)
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
    var_out = nc.createVariable(var, "f", dimensions=("time", "y", "x"), fill_value=-2e9, zlib=True, complevel=2)
    var_out.units = "Celsius"
    var_out.long_name = "theta_ocean"
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = np.repeat(theta, m * n).reshape(nt, n, m)

    var = "salinity_ocean"
    var_out = nc.createVariable(var, "f", dimensions=("time", "y", "x"), fill_value=-2e9, zlib=True, complevel=2)
    var_out.units = "g/kg"
    var_out.long_name = "salinity_ocean"
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = np.repeat(salinity, m * n).reshape(nt, n, m)

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


col_dict = {
    "ICES": "#6baed6",
    "GINR": "#c6dbef",
    "GINR S26": "#006d2c",
    "OMG Fjord": "#54278f",
    "OMG Bay": "#08519c",
    "XCTD Fjord": "#9e9ac8",
    "Fjord": "#9e9ac8",
    "Bay": "#6baed6",
}
ms = 2
mew = 0.25

fontsize = 6
lw = 0.65
aspect_ratio = 0.35
markersize = 2

params = {
    "backend": "ps",
    "axes.linewidth": 0.25,
    "lines.linewidth": lw,
    "axes.labelsize": fontsize,
    "font.size": fontsize,
    "xtick.direction": "in",
    "xtick.labelsize": fontsize,
    "xtick.major.size": 2.5,
    "xtick.major.width": 0.25,
    "ytick.direction": "in",
    "ytick.labelsize": fontsize,
    "ytick.major.size": 2.5,
    "ytick.major.width": 0.25,
    "legend.fontsize": fontsize,
    "font.size": fontsize,
}

plt.rcParams.update(params)

if __name__ == "__main__":

    # depths to average over
    depth_min = 225
    depth_max = 275
    # depth for freezing point calculation
    depth = 250
    salinity = 34
    grid_spacing = 18000

    training_iterations = 1000

    # Choose the temporal averaging window. Using "1W" instead of "1D" produces much smoother results
    freq = "1D"

    start_date = datetime(1980, 1, 1)
    end_date = datetime(2021, 1, 1)
    end_date_yearly = datetime(2021, 1, 2)

    calendar = "standard"
    units = "days since 1980-1-1"
    cdftime_days = utime(units, calendar)

    # create list with dates from start_date until end_date with
    # periodicity prule.
    bnds_datelist = list(rrule.rrule(rrule.DAILY, dtstart=start_date, until=end_date_yearly))
    bnds_datelist = list(rrule.rrule(rrule.WEEKLY, dtstart=start_date, until=end_date_yearly))
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

    # Create monthly forcing
    dates = pd.date_range(start=start_date, end=end_date, freq="1W")
    X_new = [to_decimal_year(date) for date in dates]
    X_test = torch.tensor(X_new).to(torch.float)

    # We could use "init" to create tie points, in particular, force to a point in 1980-1-1
    init = pd.read_csv("init/init.csv", parse_dates=["Date"])

    ginr = pd.read_csv("ginr/ginr_disko_bay_250m.csv", parse_dates=["Date"])
    ginr = ginr.set_index("Date").drop(columns=["Unnamed: 0"])
    ginr = ginr.groupby(pd.Grouper(freq=freq)).mean().dropna(subset=["Temperature [Celsius]", "Salinity [g/kg]"])

    ginr_ctd26 = pd.read_csv("ginr/ginr_ctd_station_26.csv").dropna()
    ginr_s26_S = pd.read_csv("ginr/GINR-S26-Salinity.csv", names=["Year", "Salinity [g/kg]"])
    ginr_s26_T = pd.read_csv("ginr/GINR-S26-Temperature.csv", names=["Year", "Temperature [Celsius]"])

    omg_fjord = pd.read_csv("omg/omg_axctd_ilulissat_fjord_10s_mean_250m.csv", parse_dates=["Date"])
    omg_fjord = omg_fjord.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_fjord = (
        omg_fjord.groupby(pd.Grouper(freq=freq)).mean().dropna(subset=["Temperature [Celsius]", "Salinity [g/kg]"])
    )

    omg_bay = pd.read_csv("omg/omg_axctd_disko_bay_10s_mean_250m.csv", parse_dates=["Date"])
    omg_bay = omg_bay.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_bay = omg_bay.groupby(pd.Grouper(freq=freq)).mean().dropna(subset=["Temperature [Celsius]", "Salinity [g/kg]"])

    ices = pd.read_csv("ices/ices_disko_bay_250m.csv", parse_dates=["Date"])
    ices = ices.set_index("Date")
    ices = ices.groupby(pd.Grouper(freq=freq)).mean().dropna(subset=["Temperature [Celsius]", "Salinity [g/kg]"])

    xctd_fjord = pd.read_csv("xctd_fjord/xctd_ilulissat_fjord_250m.csv", parse_dates=["Date"])
    xctd_fjord = xctd_fjord.set_index("Date")
    xctd_fjord = (
        xctd_fjord.groupby(pd.Grouper(freq=freq)).mean().dropna(subset=["Temperature [Celsius]", "Salinity [g/kg]"])
    )

    xctd_bay = pd.read_csv("moorings/xctd_mooring_disko_bay.csv", parse_dates=["Date"])
    xctd_bay = xctd_bay.set_index("Date")
    xctd_bay = (
        xctd_bay.groupby(pd.Grouper(freq=freq)).mean().dropna(subset=["Temperature [Celsius]", "Salinity [g/kg]"])
    )

    X_init = init["Year"].values.reshape(-1, 1)
    X_ginr = ginr["Year"].values.reshape(-1, 1)
    X_ginr_ctd26 = ginr_ctd26["Year"].values.reshape(-1, 1)
    X_ginr_s26_S = ginr_s26_S["Year"].values.reshape(-1, 1)
    X_ginr_s26_T = ginr_s26_T["Year"].values.reshape(-1, 1)
    X_ices = ices["Year"].values.reshape(-1, 1)
    X_omg_bay = omg_bay["Year"].values.reshape(-1, 1)
    X_omg_fjord = omg_fjord["Year"].values.reshape(-1, 1)
    X_xctd_bay = xctd_bay["Year"].values.reshape(-1, 1)
    X_xctd_fjord = xctd_fjord["Year"].values.reshape(-1, 1)

    T_init = init["Temperature [Celsius]"].values
    T_ginr = ginr["Temperature [Celsius]"].values
    T_ginr_ctd26 = ginr_ctd26["Temperature [Celsius]"].values
    T_ginr_s26 = ginr_s26_T["Temperature [Celsius]"].values
    T_ices = ices["Temperature [Celsius]"].values
    T_omg_bay = omg_bay["Temperature [Celsius]"].values
    T_omg_fjord = omg_fjord["Temperature [Celsius]"].values
    T_xctd_bay = xctd_bay["Temperature [Celsius]"].values
    T_xctd_fjord = xctd_fjord["Temperature [Celsius]"].values

    S_init = init["Salinity [g/kg]"].values
    S_ginr = ginr["Salinity [g/kg]"].values
    S_ginr_s26 = ginr_s26_S["Salinity [g/kg]"].values
    S_ices = ices["Salinity [g/kg]"].values
    S_omg_bay = omg_bay["Salinity [g/kg]"].values
    S_omg_fjord = omg_fjord["Salinity [g/kg]"].values
    S_xctd_bay = xctd_bay["Salinity [g/kg]"].values
    S_xctd_fjord = xctd_fjord["Salinity [g/kg]"].values

    # Noise is just inferred from the data.  This is possible because there are multiple simultaneous
    # entries for some of the observations, and also because the different tasks are correlated.
    likelihood = gpytorch.likelihoods.GaussianLikelihood()

    all_data_ind = {
        "Temperature [Celsius]": {
            "GINR": {"X": X_ginr, "Y": T_ginr},
            "GINR S26": {"X": X_ginr_s26_T, "Y": T_ginr_s26},
            "ICES": {"X": X_ices, "Y": T_ices},
            "OMG Bay": {"X": X_omg_bay, "Y": T_omg_bay},
            "XCTD Fjord": {"X": X_xctd_fjord, "Y": T_xctd_fjord},
            "OMG Fjord": {"X": X_omg_fjord, "Y": T_omg_fjord},
        },
        "Salinity [g/kg]": {
            "GINR": {"X": X_ginr, "Y": S_ginr},
            "GINR S26": {"X": X_ginr_s26_S, "Y": S_ginr_s26},
            "ICES": {"X": X_ices, "Y": S_ices},
            "OMG Bay": {"X": X_omg_bay, "Y": S_omg_bay},
            "XCTD Fjord": {"X": X_xctd_fjord, "Y": S_xctd_fjord},
            "OMG Fjord": {"X": X_omg_fjord, "Y": S_omg_fjord},
        },
    }

    X_bay_S = np.vstack([X_ginr, X_ginr_s26_S, X_ices, X_omg_bay, X_init])
    X_bay_T = np.vstack([X_ginr, X_ginr_s26_T, X_ices, X_omg_bay, X_init])
    X_fjord = np.vstack([X_xctd_fjord, X_omg_fjord])

    T_bay = np.hstack([T_ginr, T_ginr_s26, T_ices, T_omg_bay, T_init])
    T_fjord = np.hstack([T_xctd_fjord, T_omg_fjord])
    T = np.hstack([T_bay, T_fjord])

    S_bay = np.hstack([S_ginr, S_ginr_s26, S_ices, S_omg_bay, S_init])
    S_fjord = np.hstack([S_xctd_fjord, S_omg_fjord])
    S = np.hstack([S_bay, S_fjord])

    X_bay_S_2009_2020 = X_bay_S[X_bay_S > 2009]
    X_bay_T_2009_2020 = X_bay_T[X_bay_T > 2009]
    X_fjord_2009_2020 = X_fjord[X_fjord > 2009]

    T_bay_2009_2020 = T_bay[X_bay_T.ravel() > 2009]
    T_fjord_2009_2020 = T_fjord[X_fjord.ravel() > 2009]
    T_mean_diff = T_bay_2009_2020.mean() - T_fjord_2009_2020.mean()

    S_fjord_2009_2020 = S_fjord[X_fjord.ravel() > 2009]
    S_bay_2009_2020 = S_bay[X_bay_S.ravel() > 2009]
    S_mean_diff = S_bay_2009_2020.mean() - S_fjord_2009_2020.mean()

    S_bay_n = (S_bay - S.mean()) / S.std()
    S_fjord_n = (S_fjord - S.mean()) / S.std()

    T_bay_n = (T_bay - T.mean()) / T.std()
    T_fjord_n = (T_fjord - T.mean()) / T.std()

    all_data_cat = {
        "Temperature [Celsius]": {
            "Bay": {"X": X_bay_T, "Y": T_bay},
            "Fjord": {"X": X_fjord, "Y": T_fjord},
        },
        "Salinity [g/kg]": {
            "Bay": {"X": X_bay_S, "Y": S_bay},
            "Fjord": {"X": X_fjord, "Y": S_fjord},
        },
    }

    norm_all_data_cat = {
        "Temperature [Celsius]": {
            "Bay": {"X": X_bay_T, "Y": T_bay_n},
            "Fjord": {"X": X_fjord, "Y": T_fjord_n},
        },
        "Salinity [g/kg]": {
            "Bay": {"X": X_bay_S, "Y": S_bay_n},
            "Fjord": {"X": X_fjord, "Y": S_fjord_n},
        },
    }

    norm_cat = {
        "Temperature [Celsius]": {"mean": T.mean(), "std": T.std()},
        "Salinity [g/kg]": {"mean": S.mean(), "std": S.std()},
    }

    fig, ax = plt.subplots(
        2,
        1,
        sharex="col",
        figsize=[3.2, 3.2],
        num="prognostic_all",
        clear=True,
    )
    fig.subplots_adjust(hspace=0.1)

    idx = 0
    all_samples = {}
    all_Y_pred = {}
    for key, data in all_data_cat.items():
        print(f"Training {key}")

        # Put them all together
        full_train_i = torch.cat(
            [torch.full_like(torch.tensor(data[d]["X"]), dtype=torch.long, fill_value=i) for i, d in enumerate(data)]
        )
        full_train_x = torch.cat([torch.tensor(data[d]["X"]).to(torch.float) for d in data])
        full_train_y = torch.cat([torch.tensor(data[d]["Y"]).to(torch.float) for d in data])

        # Here we have two iterms that we're passing in as train_inputs
        num_tasks = len(data)
        model = MultitaskGPModel((full_train_x, full_train_i), full_train_y, likelihood, num_tasks)

        # Find optimal model hyperparameters
        model.train()
        likelihood.train()

        # Use the adam optimizer
        optimizer = torch.optim.Adam(
            [
                {"params": model.parameters()},  # Includes GaussianLikelihood parameters
            ],
            lr=0.1,
        )

        # "Loss" for GPs - the marginal log likelihood
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

        for i in range(training_iterations):
            optimizer.zero_grad()
            output = model(full_train_x, full_train_i)
            loss = -mll(output, full_train_y)
            loss.backward()
            if (i % 50) == 0:
                print(f"  Iter %d / {training_iterations} - Loss: %.3f" % (i + 1, loss.item()))
            optimizer.step()

        # Set into eval mode
        model.eval()
        likelihood.eval()

        test_i = {d: torch.full_like(X_test, dtype=torch.long, fill_value=i) for (i, d) in enumerate(data)}

        # Make predictions---one task at a time
        # We control the task we care about using the indices

        # The gpytorch.settings.fast_pred_var flag activates LOVE (for fast variances)
        # See https://arxiv.org/abs/1803.06058

        # The sampling from the model produces:
        """
        NumericalWarning: Runtime Error when computing Cholesky decomposition: 
        Matrix not positive definite after repeatedly adding jitter up to 1.0e-04. 
        Original error on first attempt: cholesky_cpu: U(5,5) is zero, singular U.
        Using RootDecomposition.
        """
        n_samples = 10
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            print(f"{key}: calculating mean")
            Y_pred = {d: likelihood(model(X_test, test_i[d])) for d in test_i}
            all_Y_pred[key] = Y_pred
            print(f"{key}: sampling")
            samples = {
                d: model(X_test, test_i[d]).sample(
                    sample_shape=torch.Size(
                        [
                            n_samples,
                        ]
                    )
                )
                for d in test_i
            }
            all_samples[key] = samples

        for k, v in Y_pred.items():

            lower, upper = v.confidence_region()
            ax[idx].fill_between(
                X_test.numpy().squeeze(),
                lower.numpy(),
                upper.numpy(),
                color=col_dict[k],
                alpha=0.40,
                linewidth=0.25,
                label=f"{k} Sample $\pm$1-$\sigma$",
            )
            # plot the data and the true latent function
            ax[idx].plot(
                data[k]["X"],
                data[k]["Y"],
                "o",
                color=col_dict[k],
                ms=ms,
                mec="k",
                mew=mew,
                label=f"{k} Observation",
                alpha=0.75,
            )
            ax[idx].plot(
                X_test.numpy(),
                v.mean.numpy().T,
                color=col_dict[k],
                linewidth=1.0,
                label=f"{k} Sample Mean",
            )

        ax[idx].set_ylabel(key)

        idx += 1

    # Use the Fjord GP for Temperature and Bay for Salinity, but correct using the difference in the 2009-2020 mean
    for s, (temperature, salinity) in enumerate(
        zip(
            all_samples["Temperature [Celsius]"]["Fjord"].numpy(),
            all_samples["Salinity [g/kg]"]["Fjord"].numpy(),
        )
    ):
        # salinity -= S_mean_diff
        # temperature -= T_mean_diff
        if s == 0:
            ax[0].plot(X_new, temperature, color=col_dict["Fjord"], linewidth=0.2, label=f"{k} Sample")
            ax[1].plot(X_new, salinity, color=col_dict["Fjord"], linewidth=0.2, label=f"{k} Sample")
        else:
            ax[0].plot(X_new, temperature, color=col_dict["Fjord"], linewidth=0.2)
            ax[1].plot(X_new, salinity, color=col_dict["Fjord"], linewidth=0.2)
        theta_ocean = temperature - melting_point_temperature(depth, salinity)
        ofile = f"jib_ocean_forcing_{s}_1980_2020.nc"
        create_nc(ofile, theta_ocean, salinity, grid_spacing, start_date, end_date, calendar, units)

    salinity_mean = all_Y_pred["Salinity [g/kg]"]["Bay"].mean.numpy()
    temperature_mean = all_Y_pred["Temperature [Celsius]"]["Fjord"].mean.numpy()
    salinity_mean -= S_mean_diff
    temperature_mean -= 0
    # ax[0].plot(X_new, temperature_mean, color="0", linewidth=1.0)
    # ax[1].plot(X_new, salinity_mean, color=col_dict["Fjord"], linewidth=0.75, linestyle="dashed")
    theta_ocean = temperature - melting_point_temperature(depth, salinity)
    ofile = f"jib_ocean_forcing_ctrl_1980_2020.nc"
    create_nc(ofile, theta_ocean, salinity, grid_spacing, start_date, end_date, calendar, units)

    ax[1].set_xlabel("Year")
    ax[1].set_xlim(1980, 2021)
    ax[0].set_ylim(0, 5)
    ax[1].set_ylim(33, 35)
    handles, labels = ax[0].get_legend_handles_labels()
    m_handles = [handles[2], handles[3], handles[6], handles[4], handles[0], handles[1], handles[5]]
    m_labels = [labels[2], labels[3], labels[6], labels[4], labels[0], labels[1], labels[5]]
    l1 = ax[0].legend(m_handles, m_labels, ncol=2)
    l1.get_frame().set_linewidth(0.0)
    l1.get_frame().set_alpha(0.0)

    set_size(3.35, 3.35)
    fig.savefig("jib_ocean_forcing_1980_2020.pdf")

    fig, ax = plt.subplots(
        2,
        1,
        sharex="col",
        figsize=[6.2, 6.2],
        num="prognostic_all",
        clear=True,
    )
    fig.subplots_adjust(hspace=0.1)

    idx = 0
    for key, data in all_data_ind.items():

        for k, v in data.items():

            ax[idx].plot(data[k]["X"], data[k]["Y"], "o", color=col_dict[k], ms=ms, mec="k", mew=mew, label=k)
            ax[idx].set_ylabel(key)

        idx += 1

    ax[1].set_xlabel("Year")
    ax[1].set_xlim(1980, 2021)
    ax[0].set_ylim(0, 5)
    ax[1].set_ylim(33, 35)
    legend = ax[0].legend(loc="upper left", ncol=2)
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.0)

    set_size(3.35, 3.35)
    fig.savefig("jib_ocean_observation_1980_2020.pdf")
