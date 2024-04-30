#!/usr/bin/env python
# Copyright (C) 2020-21 Andy Aschwanden

from datetime import datetime
import numpy as np
import pandas as pd
import pylab as plt
import xarray as xr

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
    """Convert datetime date to decimal year"""
    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction

col_dict = {
    "ICES": "#6baed6",
    "GINR": "#deebf7",
    "OMG (Fjord)": "#005a32",
    "OMG": "#084594",
    "XCTD (Fjord)": "#74c476",
}

# col_dict = {
#     "ICES (Bay)": "#e41a1c",
#     "OMG": "#e41a1c",
#     "GINR (Bay)": "#377eb8",
#     "GINR": "#377eb8",
#     "OMG (Fjord)": "#4daf4a",
#     "OMG (Bay)": "#984ea3",
#     "ICES": "#984ea3",
#     "XCTD (Fjord)": "#ff7f00",
# }

ms = 4
mew = 0.25

fontsize = 6
lw = 0.65
aspect_ratio = 0.35

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

def melting_point_temperature(depth, salinity):
    a = [-0.0575, 0.0901, -7.61e-4]
    return a[0] * salinity + a[1] + a[2] * depth

depth = 250

plt.rcParams.update(params)

if __name__ == "__main__":


    # Choose the temporal averaging window. Using "1W" instead of "1D" produces much smoother results
    freq = "1D"

    ginr = pd.read_csv("ginr/ginr_disko_bay_250m.csv", parse_dates=["Date"])
    ginr = ginr.set_index("Date").drop(columns=["Unnamed: 0"])
    ginr = ginr.groupby(pd.Grouper(freq=freq)).mean().dropna(subset=["Temperature [Celsius]", "Salinity [g/kg]"])

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

    X_ginr = ginr["Year"].values.reshape(-1, 1)
    X_ginr_s26_T = ginr_s26_T["Year"].values.reshape(-1, 1)

    X_ginr_all = np.vstack([X_ginr, X_ginr_s26_T])
    X_ices = ices["Year"].values.reshape(-1, 1)
    X_omg_bay = omg_bay["Year"].values.reshape(-1, 1)
    X_omg_fjord = omg_fjord["Year"].values.reshape(-1, 1)
    X_xctd_bay = xctd_bay["Year"].values.reshape(-1, 1)
    X_xctd_fjord = xctd_fjord["Year"].values.reshape(-1, 1)
    X = np.vstack([X_ginr, X_ginr_s26_T, X_ices, X_omg_bay]).ravel()


    S_ginr = ginr["Salinity [g/kg]"].values
    S_ginr_s26 = ginr_s26_S["Salinity [g/kg]"].values
    S_ginr_all = np.hstack([S_ginr, S_ginr_s26])
    S_ices = ices["Salinity [g/kg]"].values
    S_omg_bay = omg_bay["Salinity [g/kg]"].values

    T_ginr = ginr["Temperature [Celsius]"].values
    T_ginr_s26 = ginr_s26_T["Temperature [Celsius]"].values
    T_ices = ices["Temperature [Celsius]"].values
    T_omg_bay = omg_bay["Temperature [Celsius]"].values
    
    T_ginr -= melting_point_temperature(depth, S_ginr) 
    T_ices -= melting_point_temperature(depth, S_ices)
    T_omg_bay -= melting_point_temperature(depth, S_omg_bay)
    
    T_omg_fjord = omg_fjord["Temperature [Celsius]"].values
    T_xctd_bay = xctd_bay["Temperature [Celsius]"].values
    T_xctd_fjord = xctd_fjord["Temperature [Celsius]"].values
    T = np.hstack([T_ginr, T_ginr_s26, T_ices, T_omg_bay])
    all_data_ind = {
        "GINR": {"X": X_ginr, "Y": T_ginr},
        "ICES": {"X": X_ices, "Y": T_ices},
        "OMG": {"X": X_omg_bay, "Y": T_omg_bay},
    }

    df = pd.concat(
        [
            pd.DataFrame(
                data=np.hstack([v["Y"].reshape(-1, 1), np.repeat(k, len(v["Y"])).reshape(-1, 1)]),
                index=v["X"].ravel(),
                columns=["Temperature [Celsius]", "Dataset"],
            )
            for k, v in all_data_ind.items()
        ]
    )
    df.to_csv("jib_observations.csv")

        
    
    fig, ax = plt.subplots(
        1,
        1,
        clear=True,
    )

    gcms = ["ACCESS1-3_rcp85", "CNRM-CM6_ssp126", "CNRM-CM6_ssp585", "CNRM-ESM2_ssp585", "CSIRO-Mk3.6_rcp85", "HadGEM2-ES_rcp85",  "IPSL-CM5-MR_rcp85", "MIROC-ESM-CHEM_rcp26", "MIROC-ESM-CHEM_rcp85", "NorESM1-M_rcp85", "UKESM1-CM6_ssp585"]
    gcms = ["ACCESS1-3_rcp85", "CNRM-CM6_ssp126", "CNRM-ESM2_ssp585", "CSIRO-Mk3.6_rcp85", "HadGEM2-ES_rcp85",  "IPSL-CM5-MR_rcp85", "MIROC-ESM-CHEM_rcp85", "NorESM1-M_rcp85", "UKESM1-CM6_ssp585"]
    gcm_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    labels = []
    handles = []
    for k, gcm in enumerate(gcms):
        print(f"Adding {gcm}")
        with xr.open_dataset(f"../ismip6/MAR3.9_{gcm}_ocean_1960-2100_v4.nc") as ds:
            time = ds["time"]
            year = [to_decimal_year(d) for d in time.to_series()]
            theta = ds.sel(x=[-187924], y=[-2272156], method="nearest").stack(dim=["x", "y"])["theta_ocean"]
            h = ax.plot(year, theta, "d", color=gcm_colors[k], ms=3, mec="k", mew=mew, alpha=0.5)
            handles.append(h[0])
            labels.append(gcm)
    for data_set, data in all_data_ind.items():
        ax.plot(
            data["X"],
            data["Y"],
            "o",
            color=col_dict[data_set],
            ms=ms,
            mec="k",
            mew=mew,
            label=f"{data_set}",
        )
    ax.set_xlabel("Year")
    ax.set_ylabel("Theta (Celsius)")
    ax.set_xlim(1980, 2021)
    ax.set_ylim(0, 10)
    l1 = ax.legend(loc="upper left", title="Observations")
    l1.get_frame().set_linewidth(0.0)
    l1.get_frame().set_alpha(0.0)
    l2 = ax.legend(handles=handles, labels=labels, loc="upper right", ncols=2, title="Simulations")
    l2.get_frame().set_linewidth(0.0)
    l2.get_frame().set_alpha(0.0)
    ax.add_artist(l1)
    ax.add_artist(l2)
    set_size(3.2, 2.2)

    fig.savefig("jib_obs_1980_2020.pdf", pad_inches="tight")
