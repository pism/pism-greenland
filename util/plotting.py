#!/usr/bin/env python

# Copyright (C) 2019 Andy Aschwanden

from argparse import ArgumentParser
from netCDF4 import Dataset as NC
from netCDF4 import num2date

import cf_units
import numpy as np
import os
import pylab as plt
import re
from pypismtools.pypismtools import smooth

from datetime import datetime

try:
    from pypismtools import unit_converter
except:
    from pypismtools.pypismtools import unit_converter


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

plot_var = "grounding_line_flux"
mass_ounits = "Gt"
flux_ounits = "Gt year-1"

vars_dict = {
    "grounding_line_flux": {
        "ounits": "Gt year-1",
        "vtype": "mass",
        "sign": 1,
        "as19_name": "grounding_line_flux",
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

fig = plt.figure()
ax = fig.add_subplot(111)

for k, ifile in enumerate(ifiles):
    nc = NC(ifile)
    exp = re.search("id_(.+?)_", ifile).group(1).upper()
    time = nc.variables["time"]
    time_units = time.units
    time_calendar = time.calendar
    date = num2date(time[:], units=time_units, calendar=time_calendar)
    var_vals = nc.variables[plot_var][:]
    iunits = nc.variables[plot_var].units
    # ax.plot_date(date, var_vals, "-", color=colors[k], linestyle="solid", linewidth=0.2)
    var_vals_smoothed = smooth(var_vals, 13 * 30)
    ax.plot_date(date, var_vals_smoothed, "-", color=colors[k], linestyle="solid", linewidth=1.0, label=exp)
    nc.close()
legend = ax.legend()
legend.get_frame().set_linewidth(0.0)
legend.get_frame().set_alpha(0.0)

# ax.set_xlim(datetime(2015, 1, 1), datetime(2100, 1, 1))
ax.set_xlabel("Year")
ax.set_ylabel(ylabel)
set_size(6, 3)
ofile = "calving_{}.pdf".format(plot_var)
print("  saving to {}".format(ofile))
fig.savefig(ofile, bbox_inches="tight")
