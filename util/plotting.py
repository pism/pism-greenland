#!/usr/bin/env python

# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
import xarray as xr

import cf_units
import cftime
import numpy as np
import os
import pandas as pd
import pylab as plt
import re
from pypismtools.pypismtools import smooth

from datetime import datetime

try:
    from pypismtools import unit_converter
except:
    from pypismtools.pypismtools import unit_converter


def set_size(w, h, ax=None):
    """w, h: width, height in inches"""

    if not ax:
        ax = plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w) / (r - l)
    figh = float(h) / (t - b)
    ax.figure.set_size_inches(figw, figh)


# Set up the option parser
parser = ArgumentParser()
parser.description = "A script for PISM output files to time series plots using pylab/matplotlib."
parser.add_argument("FILE", nargs="*")
options = parser.parse_args()
ifiles = options.FILE
# Conversion between giga tons (Gt) and millimeter sea-level equivalent (mmSLE)
gt2mmSLE = 1.0 / 362.5
gt2cmSLE = 1.0 / 362.5 / 10.0
gt2mSLE = 1.0 / 362.5 / 1000.0

plot_var = "total_grounding_line_flux"
# plot_var = "tendency_of_ice_mass_due_to_discharge"
# plot_var = "tendency_of_ice_mass_due_to_basal_mass_flux"
mass_ounits = "Gt"
flux_ounits = "Gt year-1"

vars_dict = {
    "total_grounding_line_flux": {
        "ounits": "Gt year-1",
        "vtype": "mass",
        "sign": 1,
        "as19_name": "total_grounding_line_flux",
        "normalize": False,
        "ylabel": "flux (Gt/yr)",
        "mass2sle": 1,
    },
    "grounding_line_flux": {
        "ounits": "Gt year-1",
        "vtype": "mass",
        "sign": 1,
        "as19_name": "total_grounding_line_flux",
        "normalize": False,
        "ylabel": "flux (Gt/yr)",
        "mass2sle": 1,
    },
    "limnsw": {
        "ounits": "Gt",
        "vtype": "mass",
        "sign": -1,
        "as19_name": "limnsw",
        "normalize": True,
        "ylabel": "sea level contribution (cm)",
        "mass2sle": gt2cmSLE,
    },
    "tendacabf": {
        "ounits": "Gt year-1",
        "vtype": "flux",
        "sign": 1,
        "as19_name": "tendency_of_ice_mass_due_to_surface_mass_balance",
        "normalize": False,
        "ylabel": "flux (Gt/yr)",
        "mass2sle": 1,
    },
    "tendlicalvf": {
        "ounits": "Gt year-1",
        "vtype": "flux",
        "sign": 1,
        "as19_name": "tendency_of_ice_mass_due_to_discharge",
        "normalize": False,
        "ylabel": "flux (Gt/yr)",
        "mass2sle": 1,
    },
    "tendency_of_ice_mass_due_to_discharge": {
        "ounits": "Gt year-1",
        "vtype": "flux",
        "sign": 1,
        "as19_name": "tendency_of_ice_mass_due_to_discharge",
        "normalize": False,
        "ylabel": "flux (Gt/yr)",
        "mass2sle": 1,
    },
    "tendency_of_ice_mass_due_to_basal_mass_flux": {
        "ounits": "Gt year-1",
        "vtype": "flux",
        "sign": 1,
        "as19_name": "tendency_of_ice_mass_due_to_basal_mass_flux",
        "normalize": False,
        "ylabel": "flux (Gt/yr)",
        "mass2sle": 1,
    },
    "tendency_of_ice_mass_due_to_calving": {
        "ounits": "Gt year-1",
        "vtype": "flux",
        "sign": 1,
        "as19_name": "tendency_of_ice_mass_due_to_calving",
        "normalize": False,
        "ylabel": "flux (Gt/yr)",
        "mass2sle": 1,
    },
}

rcp_col_dict = {"CTRL": "k", "85": "#990002", "45": "#5492CD", "26": "#003466"}
rcp_gray_col_dict = {"CTRL": "k", "85": "0.70", "45": "0.75", "26": "0.80"}
rcp_shade_col_dict = {"CTRL": "k", "85": "#F4A582", "45": "#92C5DE", "26": "#4393C3"}
rcp_dict = {"26": "RCP 2.6", "45": "RCP 4.5", "85": "RCP 8.5", "CTRL": "CTRL"}

colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

print("Plotting {}".format(plot_var))

m_var = vars_dict[plot_var]
sign = m_var["sign"]
as19_name = m_var["as19_name"]
normalize = m_var["normalize"]
ylabel = m_var["ylabel"]
mass2sle = m_var["mass2sle"]

D = pd.read_csv("~/Google Drive/My Drive/data/mankoff_discharge/gate_merged.csv", parse_dates=["Date"])
D = D[D["Gate"] == 184]
fig = plt.figure()
ax = fig.add_subplot(111)

ax.errorbar(D["Date"], -D["Discharge [Gt/yr]"], yerr=D["Discharge Error [Gt/yr]"], c="#238b45", lw=1.0, elinewidth=0.5)
for k, ifile in enumerate(ifiles):
    nc = xr.open_dataset(ifile)
    print(ifile)
    exp = re.search("id_(.+?)_", ifile).group(1).upper()
    time = nc.variables["time"]
    # time_units = time.units
    # time_calendar = time.calendar
    # date = cftime.num2pydate(time, time_units)
    date = time
    data = nc.variables[plot_var]
    var_vals = np.squeeze(data)
    iunits = data.attrs["units"]
    # ax.plot_date(date, var_vals, "-", color=colors[k], linestyle="solid", linewidth=0.2)
    var_vals_smoothed = smooth(var_vals, 13)
    #    ax.plot_date(date, var_vals_smoothed, "-", linewidth=1.0, label=exp)
    ax.plot_date(date, var_vals_smoothed, "-", color="0.5", linewidth=0.5, label=exp)
    nc.close()
# legend = ax.legend()
# legend.get_frame().set_linewidth(0.0)
# legend.get_frame().set_alpha(0.0)

# ax.set_xlim(datetime(2015, 1, 1), datetime(2100, 1, 1))
ax.set_xlabel("Year")
ax.set_ylabel(ylabel)
set_size(6, 3)
ofile = "calving_{}.pdf".format(plot_var)
print("  saving to {}".format(ofile))
fig.savefig(ofile, bbox_inches="tight")


# df_meta = pd.read_csv("/Users/andy/Downloads/gate_meta.csv")
# df_D = pd.read_csv("/Users/andy/Downloads/gate_D.csv")
# df_D.index = df_D["Date"]
# df_D = df_D[df_D.index.isin(dates_merged)]
# df_D_nd = df_D.drop(columns=["Date"])
# df_err = pd.read_csv("/Users/andy/Downloads/gate_err.csv")
# df_err.index = df_err["Date"]
# df_err = df_err[df_err.index.isin(dates_merged)]
# df_err_nd = df_err.drop(columns=["Date"])
# dfs = [pd.DataFrame(data = np.hstack([df_D["Date"].values.reshape(-1,1), np.repeat(col, len(df_D_nd)).reshape(-1, 1), df_D_nd[col].values.reshape(-1,1), df_err_nd[col].values.reshape(-1,1)]), columns=["Date", "Gate", "Discharge [Gt/yr]", "Discharge Error [Gt/yr]"]).astype({"Gate": int}) for col in df_D_nd]
# all_D = pd.concat(dfs)
# df = pd.merge(all_D, df_meta, left_on="Gate", right_on="gate").drop(columns=["gate"]).astype({"Date": "datetime64[ns]"})
# df.to_csv("/Users/andy/Downloads/gate_merged.csv")
