#!/usr/bin/env python

# Copyright (C) 2019-21 Andy Aschwanden

from argparse import ArgumentParser
from netCDF4 import Dataset as NC

import nc_time_axis
import cf_units
import cftime
import numpy as np
import os
import pylab as plt
from glob import glob
from pathlib import Path
import re

try:
    from pypismtools import unit_converter, trend_estimator
except:
    from pypismtools.pypismtools import unit_converter, trend_estimator


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


def smooth(x, window_len=11, window="hanning"):
    """
    Smooth the data using a window with requested size (running mean,
    moving average, low pass filtering).

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x : array_like, the input signal
    window_len : the dimension of the smoothing window; should be an odd integer
    window : the type of window from "flat", "hanning", "hamming",
    "bartlett", "blackman" flat window will produce a moving average smoothing.

    Returns
    -------
    y : the smoothed signal

    Example
    -------
    t = np.linspace(-2,2,0.1)
    x = np.sin(t) + np.randn(len(t))*0.1
    y = smooth(x)

    See also
    --------
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve
    scipy.signal.lfilter

    Notes
    -----
    Downloaded from http://www.scipy.org/Cookbook/SignalSmooth.

    TODO
    ----
    the window parameter could be the window itself if an array instead
    of a string
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )

    s = np.r_[2 * x[0] - x[window_len:1:-1], x, 2 * x[-1] - x[-1:-window_len:-1]]

    if window == "flat":  # moving average
        w = np.ones(window_len, "d")
    else:
        w = eval("np." + window + "(window_len)")

    y = np.convolve(w / w.sum(), s, mode="same")
    return y[window_len - 1 : -window_len + 1]


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

nc0 = NC(ifiles[0])
profiles = nc0.variables["profile_name"][:]
v_units = nc0.variables[var].units
nc0.close()

for p_id, profile in enumerate(profiles):
    print(profile)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for m_file in ifiles:
        nc = NC(m_file)
        exp = re.search("id_(.+?)_", m_file).group(1).upper()
        time = nc.variables["time"]
        time_units = time.units
        dates = cftime.num2pydate(time, time_units)
        v = nc.variables[var]
        data = v[p_id, :, p_id]
        # ax.plot_date(dates, data, "-", c="0.5", ms=0, lw=0.2)
        lw = 0.5
        ax.plot_date(dates, smooth(data, 13), "-", ms=0, lw=lw, label=exp)
        nc.close()

    if obs_file is not None:
        nc0 = NC(obs_file)
        profiles = nc0.variables["profile_name"][:]
        v_units = nc0.variables[var].units
        time = nc.variables["time"]
        time_units = time.units
        dates = cftime.num2pydate(time, time_units)
        v = nc.variables[var]
        data = v[p_id, :, p_id]
        lw = 1.5
        ax.plot_date(
            dates, smooth(data, 13), "o", color="k", ms=1, lw=lw, label="ITS_LIVE"
        )
        nc0.close()
        ax.set_ylim(0, data.max())

    ax.set_xlim(cftime.datetime(1980, 1, 1), cftime.datetime(2020, 1, 1))
    ax.set_ylabel(f"{var} ({v_units})")
    ax.set_title(profile)
    legend_1 = ax.legend(loc="upper left", ncol=3)
    legend_1.get_frame().set_linewidth(0.0)
    legend_1.get_frame().set_alpha(0.0)

    set_size(6, 2)
    fig.savefig(f"pt_{profile}_test.pdf")
