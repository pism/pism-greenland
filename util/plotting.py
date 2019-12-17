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

plot_var = "limnsw"
mass_ounits = "Gt"
flux_ounits = "Gt year-1"

vars_dict = {
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
}

rcp_col_dict = {"CTRL": "k", "85": "#990002", "45": "#5492CD", "26": "#003466"}
rcp_gray_col_dict = {"CTRL": "k", "85": "0.70", "45": "0.75", "26": "0.80"}
rcp_shade_col_dict = {"CTRL": "k", "85": "#F4A582", "45": "#92C5DE", "26": "#4393C3"}
rcp_dict = {"26": "RCP 2.6", "45": "RCP 4.5", "85": "RCP 8.5", "CTRL": "CTRL"}

as19_dir = "/Volumes/zachariae/pism-gris/stability/2018_09_les/scalar_ensstat/"

for plot_var in ("limnsw", "tendacabf"):

    print("Plotting {}".format(plot_var))

    m_var = vars_dict[plot_var]
    sign = m_var["sign"]
    as19_name = m_var["as19_name"]
    normalize = m_var["normalize"]
    ylabel = m_var["ylabel"]
    mass2sle = m_var["mass2sle"]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for rcp in ["85", "26"]:
        as19_vals = []
        for pctl in [5, 50, 95]:
            ifile = os.path.join(
                as19_dir, "enspctl{pctl}_gris_g1800m_v3a_rcp_{rcp}_0_1000.nc".format(pctl=pctl, rcp=rcp)
            )
            nc = NC(ifile)
            time = nc.variables["time"]
            time_units = time.units
            time_calendar = time.calendar
            date = num2date(time[7:94], units=time_units, calendar=time_calendar)
            var_vals = nc.variables[as19_name][7:94]
            if normalize:
                var_vals -= nc.variables[as19_name][7]
            iunits = nc.variables[as19_name].units
            as19_vals.append(sign * unit_converter(var_vals, iunits, m_var["ounits"]) * mass2sle)
            nc.close()
        ax.fill_between(
            date, as19_vals[0], as19_vals[2], color=rcp_gray_col_dict[rcp], linewidth=0, label="AS19 RCP{}".format(rcp)
        )
        ax.plot_date(date, as19_vals[0], "-", color="0.25", linewidth=0.2)
        ax.plot_date(date, as19_vals[2], "-", color="0.25", linewidth=0.2)
        plt.title(plot_var)

    for ifile in ifiles:
        nc = NC(ifile)
        exp = re.search("EXP-(.+?)_", ifile).group(1).upper()
        time = nc.variables["time"]
        time_units = time.units
        time_calendar = time.calendar
        date = num2date(time[:], units=time_units, calendar=time_calendar)
        var_vals = nc.variables[plot_var][:]
        if normalize:
            var_vals -= nc.variables[plot_var][0]
        iunits = nc.variables[plot_var].units
        var_vals = sign * unit_converter(var_vals, iunits, m_var["ounits"]) * mass2sle
        ax.plot_date(date, var_vals, "-", linestyle="-", linewidth=0.4, label=exp)
        nc.close()
    legend = ax.legend()
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.0)

    ax.set_xlim(datetime(2015, 1, 1), datetime(2100, 1, 1))
    ax.set_xlabel("Year")
    ax.set_ylabel(ylabel)
    set_size(6, 3)
    ofile = "UAF_PISM_{}.pdf".format(plot_var)
    print("  saving to {}".format(ofile))
    fig.savefig(ofile, bbox_inches="tight")
