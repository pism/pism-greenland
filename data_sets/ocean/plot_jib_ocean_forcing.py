#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

from cftime import utime
from dateutil import rrule
from dateutil.parser import parse
from datetime import datetime
import numpy as np
import pandas as pd
import pylab as plt
from itertools import product
import os
import time
import glob
from netCDF4 import Dataset as NC


def toDecimalYear(date):
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year + 1, month=1, day=1)
    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed / yearDuration

    return date.year + fraction


def melting_point_temperature(depth, salinity):
    a = [-0.0575, 0.0901, -7.61e-4]
    return a[0] * salinity + a[1] + a[2] * depth


if __name__ == "__main__":

    # depths to average over
    depth_min = 225
    depth_max = 275
    # depth for freezing point calculation
    depth = 250
    salinity = 34
    freq = "1D"
    xi = 526
    yi = 1176
    T_m = melting_point_temperature(depth, salinity)
    X_mar = np.arange(1960, 2101)

    ginr_df = pd.read_csv("ginr/ginr_ctd_station_26.csv")

    omg_df = pd.read_csv("omg/omg_axctd_ilulissat_fjord_10s_mean.csv", parse_dates=["Date"])
    omg_df = omg_df[(omg_df["Depth [meters]"] <= depth_max) & (omg_df["Depth [meters]"] >= depth_min)]
    omg_df = omg_df.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_df = omg_df.groupby(pd.Grouper(freq=freq)).mean().dropna()

    ices_df = pd.read_csv("ices/ices_disko_bay.csv")
    ices_df = ices_df[(ices_df["Depth [m]"] >= depth_min) & (ices_df["Depth [m]"] <= depth_max)].reset_index(drop=True)
    ices_time = pd.to_datetime(ices_df.Date, format="%Y/%m/%d %H:%M:%S")
    ices_df["Date"] = ices_time
    ices_df = ices_df.set_index("Date")
    ices_df = ices_df.groupby(pd.Grouper(freq=freq)).mean().dropna()
    ices_time = [toDecimalYear(d) for d in ices_df.index]
    ices_df["Year"] = ices_time

    xctd_if_df = pd.read_csv("xctd_fjord/xctd_ilulissat_fjord.csv", parse_dates=["Date"]).dropna()
    xctd_if_df = xctd_if_df[
        (xctd_if_df["Depth [m]"] >= depth_min) & (xctd_if_df["Depth [m]"] <= depth_max)
    ].reset_index(drop=True)
    xctd_if_df = xctd_if_df.set_index("Date")
    xctd_if_df = xctd_if_df.groupby(pd.Grouper(freq=freq)).mean().dropna()
    xctd_if_time = [toDecimalYear(d) for d in xctd_if_df.index]
    xctd_if_df["Year"] = xctd_if_time

    xctd_db_df = pd.read_csv("moorings/xctd_mooring_disko_bay.csv", parse_dates=["Date"])
    xctd_db_df = xctd_db_df.set_index("Date")
    xctd_db_df = xctd_db_df.groupby(pd.Grouper(freq=freq)).mean().dropna()
    xctd_db_time = [toDecimalYear(d) for d in xctd_db_df.index]
    xctd_db_df["Year"] = xctd_db_time

    X_ginr = ginr_df["Year"].values.reshape(-1, 1)
    X_ices = ices_df["Year"].values.reshape(-1, 1)
    X_omg = omg_df["Year"].values.reshape(-1, 1)
    X_xctd_if = xctd_if_df["Year"].values.reshape(-1, 1)
    X_xctd_db = xctd_db_df["Year"].values.reshape(-1, 1)

    y_ginr = ginr_df["Temperature [Celsius]"].values
    y_ices = ices_df["Temperature [Celsius]"].values
    y_omg = omg_df["Temperature [Celsius]"].values
    y_xctd_if = xctd_if_df["Temperature [Celsius]"].values
    y_xctd_db = xctd_db_df["Temperature [Celsius]"].values

    fig = plt.figure()
    ax = fig.gca()

    ax.plot(X_ices, y_ices, "o", color="#08519c", ms=4, mec="k", mew=0.1, label="ICES (Disko Bay)")
    ax.plot(X_ginr, y_ginr, "o", color="#6baed6", ms=4, mec="k", mew=0.1, label="GINR Station 26 (Disko Bay)")
    ax.plot(X_xctd_db, y_xctd_db, "o", color="#c6dbef", mec="k", mew=0.1, ms=4, label="Mooring (Disko Bay)")
    ax.plot(X_xctd_if, y_xctd_if, "o", color="#a50f15", ms=4, mec="k", mew=0.1, label="XCTD (Ilulissat Fjord)")
    ax.plot(X_omg, y_omg, "o", color="#fb6a4a", ms=4, mec="k", mew=0.1, label="OMG (Ilulissat Fjord)")

    ax.set_xlabel("Year")
    ax.set_ylabel("Temperature (Celsius)")
    ax.set_xlim(1980, 2021)
    ax.set_ylim(0, 4)
    legend = ax.legend()
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.0)

    fig.savefig("disko-bay-temperature-obs.pdf")
    plt.cla()
    plt.close(fig)

    fig = plt.figure()
    ax = fig.gca()

    ax.plot(X_ices, y_ices, "o", color="#08519c", ms=4, mec="k", mew=0.1, label="ICES (Disko Bay)")
    ax.plot(X_ginr, y_ginr, "o", color="#6baed6", ms=4, mec="k", mew=0.1, label="GINR Station 26 (Disko Bay)")
    ax.plot(X_xctd_db, y_xctd_db, "o", color="#c6dbef", mec="k", mew=0.1, ms=4, label="Mooring (Disko Bay)")
    ax.plot(X_xctd_if, y_xctd_if, "o", color="#a50f15", ms=4, mec="k", mew=0.1, label="XCTD (Ilulissat Fjord)")
    ax.plot(X_omg, y_omg, "o", color="#fb6a4a", ms=4, mec="k", mew=0.1, label="OMG (Ilulissat Fjord)")

    for k, m_file in enumerate(glob.glob("../ismip6/MAR3.9_*ocean_1960-2100_v4.nc")):
        if m_file.find("ctrl") == -1:
            nc = NC(m_file)
            theta = nc.variables["theta_ocean"][:, yi, xi]
            T_mar = theta + T_m
            if k == 0:
                l = ax.plot(X_mar, T_mar, color="0.5", linewidth=0.25, label="MAR")
            else:
                ax.plot(X_mar, T_mar, color="0.5", linewidth=0.25)

    ax.set_xlabel("Year")
    ax.set_ylabel("Temperature (Celsius)")
    ax.set_xlim(1980, 2021)
    ax.set_ylim(0, 4)
    legend = ax.legend()
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.0)
    fig.savefig("disko-bay-temperature-obs-mar.pdf")

    plt.cla()
    plt.close(fig)
