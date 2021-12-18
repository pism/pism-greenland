#!/usr/bin/env python

# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
import xarray as xr

import cf_units
import numpy as np
import os
import pandas as pd
import pylab as plt

from datetime import datetime


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
parser.description = (
    "A script for PISM output files to time series plots using pylab/matplotlib."
)
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

D = pd.read_csv(
    "~/Google Drive/My Drive/data/mankoff_discharge/gate_merged.csv",
    parse_dates=["Date"],
)
D = D[D["Gate"] == 184]

df = pd.read_csv("fldsum_test.csv", parse_dates=["time"])
fig = plt.figure()
ax = fig.add_subplot(111)


def plot_ts(f, ax):
    ax.plot(
        f.index,
        f["total_grounding_line_flux (Gt year-1)"],
        color="0.5",
        lw=0.25,
        alpha=0.2,
    )


[plot_ts(f, ax) for f in df.groupby(by="id").rolling(13, on="time")]

ax.fill_between(
    D["Date"],
    -D["Discharge [Gt/yr]"] - D["Discharge Error [Gt/yr]"],
    -D["Discharge [Gt/yr]"] + D["Discharge Error [Gt/yr]"],
    color="#238b45",
    alpha=0.5,
    lw=1.0,
)

ax.plot(
    D["Date"],
    -D["Discharge [Gt/yr]"],
    color="#238b45",
    lw=1.0,
    label="Mankoff",
)

ax.set_xlim(datetime(1985, 1, 1), datetime(1990, 1, 1))
ax.set_xlabel("Year")
ax.set_ylabel(ylabel)
ax.set_ylim(-75, -10)
set_size(6, 3)
ofile = "calving_{}.pdf".format(plot_var)
print("  saving to {}".format(ofile))
fig.savefig(ofile, bbox_inches="tight")
