#!/usr/bin/env python

# Copyright (C) 2019-21 Andy Aschwanden

from argparse import ArgumentParser
import xarray as xr
import numpy as np
from datetime import datetime
import os
import pylab as plt
from glob import glob
from pathlib import Path
import re


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
parser.add_argument("--obs_file", default=None)

options = parser.parse_args()
ifiles = options.FILE
obs_file = options.obs_file

var = "velsurf_mag"

nc0 = xr.open_dataset(ifiles[0])
profiles = nc0.variables["profile_name"][:]
v_units = nc0.variables[var].attrs["units"]
nc0.close()

for p_id, profile in enumerate(profiles):
    print(profile.values)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for m_file in ifiles:
        nc = xr.open_dataset(m_file)
        exp = re.search("id_(.+?)_", m_file).group(1).upper()
        dates = nc.variables["time"]
        v = nc.variables[var]
        data = v[p_id, :, p_id]
        print(m_file, data[-1].values)
        lw = 0.5
        ax.plot(dates, data, "-", color="0.5", ms=0, lw=lw)
        nc.close()

    if obs_file is not None:
        nc_obs = xr.open_dataset(obs_file)
        profiles = nc_obs.variables["profile_name"][:]
        v_units = nc_obs.variables[var].attrs["units"]
        dates = nc_obs.variables["time"]
        v_obs = nc_obs.variables[var]
        speed_obs = v_obs[p_id, :, p_id]
        v_obs_err = nc_obs.variables["v_err"]
        speed_obs_err = v_obs_err[p_id, :, p_id]
        lw = 1.5
        ax.fill_between(
            dates,
            speed_obs - 3 * speed_obs_err,
            speed_obs + 3 * speed_obs_err,
            color="#cbc9e2",
            alpha=0.5,
            lw=0.0,
        )

        ax.plot(dates, speed_obs, "o", color="#6a51a3", ms=5, lw=lw, label="ITS_LIVE")
        nc_obs.close()
        ax.set_ylim(0, 20000)

    ax.set_xlim(datetime(1985, 1, 1), datetime(2020, 1, 1))
    ax.set_ylabel(f"{var} ({v_units})")
    ax.set_title(profile.values)
    # legend_1 = ax.legend(loc="upper left", ncol=3)
    # legend_1.get_frame().set_linewidth(0.0)
    # legend_1.get_frame().set_alpha(0.0)

    set_size(6, 2)
    fig.savefig(f"pt_{profile.values}.pdf")
