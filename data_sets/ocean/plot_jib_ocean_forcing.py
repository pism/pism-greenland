#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

from datetime import datetime
import numpy as np
import pandas as pd
import pylab as plt
import glob
from netCDF4 import Dataset as NC


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

    ginr = pd.read_csv("ginr/ginr_disko_bay_250m.csv", parse_dates=["Date"])
    ginr = ginr.set_index("Date").drop(columns=["Unnamed: 0"])
    ginr = ginr.groupby(pd.Grouper(freq=freq)).mean()

    ginr_ctd26 = pd.read_csv("ginr/ginr_ctd_station_26.csv")

    omg_fjord = pd.read_csv("omg/omg_axctd_ilulissat_fjord_10s_mean_250m.csv", parse_dates=["Date"])
    omg_fjord = omg_fjord.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_fjord = omg_fjord.groupby(pd.Grouper(freq=freq)).mean()

    omg_bay = pd.read_csv("omg/omg_axctd_disko_bay_10s_mean_250m.csv", parse_dates=["Date"])
    omg_bay = omg_bay.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_bay = omg_bay.groupby(pd.Grouper(freq=freq)).mean()

    ices = pd.read_csv("ices/ices_disko_bay_250m.csv", parse_dates=["Date"])
    ices = ices.set_index("Date")
    ices = ices.groupby(pd.Grouper(freq=freq)).mean()

    xctd_fjord = pd.read_csv("xctd_fjord/xctd_ilulissat_fjord_250m.csv", parse_dates=["Date"])
    xctd_fjord = xctd_fjord.set_index("Date")
    xctd_fjord = xctd_fjord.groupby(pd.Grouper(freq=freq)).mean()

    xctd_bay = pd.read_csv("moorings/xctd_mooring_disko_bay.csv", parse_dates=["Date"])
    xctd_bay = xctd_bay.set_index("Date")
    xctd_bay = xctd_bay.groupby(pd.Grouper(freq=freq)).mean()

    X_ginr = ginr["Year"].values.reshape(-1, 1)
    X_ginr_ctd26 = ginr_ctd26["Year"].values.reshape(-1, 1)
    X_ices = ices["Year"].values.reshape(-1, 1)
    X_omg_bay = omg_bay["Year"].values.reshape(-1, 1)
    X_omg_fjord = omg_fjord["Year"].values.reshape(-1, 1)
    X_xctd_bay = xctd_bay["Year"].values.reshape(-1, 1)
    X_xctd_fjord = xctd_fjord["Year"].values.reshape(-1, 1)

    T_ginr = ginr["Temperature [Celsius]"].values
    T_ginr_ctd26 = ginr_ctd26["Temperature [Celsius]"].values
    T_ices = ices["Temperature [Celsius]"].values
    T_omg_bay = omg_bay["Temperature [Celsius]"].values
    T_omg_fjord = omg_fjord["Temperature [Celsius]"].values
    T_xctd_bay = xctd_bay["Temperature [Celsius]"].values
    T_xctd_fjord = xctd_fjord["Temperature [Celsius]"].values

    S_ginr = ginr["Salinity [g/kg]"].values
    S_ices = ices["Salinity [g/kg]"].values
    S_omg_bay = omg_bay["Salinity [g/kg]"].values
    S_omg_fjord = omg_fjord["Salinity [g/kg]"].values
    S_xctd_bay = xctd_bay["Salinity [g/kg]"].values
    S_xctd_fjord = xctd_fjord["Salinity [g/kg]"].values

    omg_bay_col = "#08519c"
    ices_bay_col = "#6baed6"
    ginr_bay_col = "#c6dbef"
    ginr_ctd26_col = "#74c476"

    omg_fjord_col = "#54278f"
    xctd_fjord_col = "#9e9ac8"
    ms = 5
    mew = 0.25
    fig, ax = plt.subplots(
        2,
        1,
        sharex="col",
        figsize=[6.2, 6.2],
        num="prognostic_all",
        clear=True,
    )
    fig.subplots_adjust(hspace=0.1)

    ax[0].plot(X_omg_bay, T_omg_bay, "o", color=omg_bay_col, ms=ms, mec="k", mew=mew, label="OMG (Disko Bay)")
    ax[0].plot(X_ices, T_ices, "o", color=ices_bay_col, ms=ms, mec="k", mew=mew, label="ICES (Disko Bay)")
    ax[0].plot(X_ginr, T_ginr, "o", color=ginr_bay_col, ms=ms, mec="k", mew=mew, label="GINR (Disko Bay)")
    ax[0].plot(
        X_ginr_ctd26, T_ginr_ctd26, "o", color=ginr_ctd26_col, ms=ms, mec="k", mew=mew, label="GINR (Station 26)"
    )
    # ax[0].plot(X_xctd_bay, T_xctd_bay, "o", color="#006d2c", mec="k", mew=mew, ms=ms, label="Mooring (Disko Bay)")
    ax[0].plot(
        X_xctd_fjord, T_xctd_fjord, "o", color=xctd_fjord_col, ms=ms, mec="k", mew=mew, label="XCTD (Ilulissat Fjord)"
    )
    ax[0].plot(
        X_omg_fjord, T_omg_fjord, "o", color=omg_fjord_col, ms=ms, mec="k", mew=mew, label="OMG (Ilulissat Fjord)"
    )

    ax[1].plot(X_omg_bay, S_omg_bay, "o", color=omg_fjord_col, ms=ms, mec="k", mew=mew, label="OMG Fjord (Disko Bay)")
    ax[1].plot(X_ices, S_ices, "o", color=ices_bay_col, ms=ms, mec="k", mew=mew, label="ICES (Disko Bay)")
    ax[1].plot(X_ginr, S_ginr, "o", color=ginr_bay_col, ms=ms, mec="k", mew=mew, label="GINR (Disko Bay)")
    # ax[1].plot(X_xctd_bay, S_xctd_bay, "o", color="#006d2c", mec="k", mew=mew, ms=ms, label="Mooring (Disko Bay)")
    ax[1].plot(
        X_xctd_fjord, S_xctd_fjord, "o", color=xctd_fjord_col, ms=ms, mec="k", mew=mew, label="XCTD (Ilulissat Fjord)"
    )
    ax[1].plot(
        X_omg_fjord, S_omg_fjord, "o", color=omg_fjord_col, ms=ms, mec="k", mew=mew, label="OMG (Ilulissat Fjord)"
    )

    ax[0].set_ylabel("Temperature (Celsius)")
    ax[0].set_ylim(0, 5)
    ax[1].set_xlabel("Year")
    ax[1].set_ylabel("Salinity (g/kg)")
    ax[1].set_xlim(1980, 2021)

    legend = ax[0].legend()
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.0)

    fig.savefig("ji_ocean_obs_all.pdf")

    for k, m_file in enumerate(glob.glob("../ismip6/MAR3.9_*ocean_1960-2100_v4.nc")):
        if m_file.find("ctrl") == -1:
            nc = NC(m_file)
            theta = nc.variables["theta_ocean"][:, yi, xi]
            T_mar = theta + T_m
            if k == 0:
                ax[0].plot(X_mar, T_mar, color="0.5", linewidth=0.25, label="MAR")
            else:
                ax[0].plot(X_mar, T_mar, color="0.5", linewidth=0.25)

    fig.savefig("disko-bay-temperature-obs-mar.pdf")

    plt.cla()
    plt.close(fig)
