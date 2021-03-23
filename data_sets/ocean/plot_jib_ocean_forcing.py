#!/usr/bin/env python
# Copyright (C) 2020-21 Andy Aschwanden

from datetime import datetime
import numpy as np
import pandas as pd
import pylab as plt


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

plt.rcParams.update(params)

if __name__ == "__main__":

    ji_mb = pd.read_csv("mass_loss/rates_JI.txt")

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

    T_ginr = ginr["Temperature [Celsius]"].values
    T_ginr_s26 = ginr_s26_T["Temperature [Celsius]"].values
    T_ginr_all = np.hstack([T_ginr, T_ginr_s26])
    T_ices = ices["Temperature [Celsius]"].values
    T_omg_bay = omg_bay["Temperature [Celsius]"].values
    T_omg_fjord = omg_fjord["Temperature [Celsius]"].values
    T_xctd_bay = xctd_bay["Temperature [Celsius]"].values
    T_xctd_fjord = xctd_fjord["Temperature [Celsius]"].values
    T = np.hstack([T_ginr, T_ginr_s26, T_ices, T_omg_bay])
    print(T.shape)
    all_data_ind = {
        "GINR": {"X": X_ginr_all, "Y": T_ginr_all},
        "ICES": {"X": X_ices, "Y": T_ices},
        "OMG": {"X": X_omg_bay, "Y": T_omg_bay},
        # "XCTD (Fjord)": {"X": X_xctd_fjord, "Y": T_xctd_fjord},
        # "OMG (Fjord)": {"X": X_omg_fjord, "Y": T_omg_fjord},
    }

    X_1980_1996 = X[X < 1997]
    X_1997_2015 = X[X >= 1997]
    T_1980_1996 = T[X < 1997]
    T_1997_2015 = T[X >= 1997]

    T_1980_1996_mean = np.mean(T_1980_1996)
    T_1997_2015_mean = np.mean(T_1997_2015)

    print(T_1997_2015_mean - T_1980_1996_mean)

    fig, ax = plt.subplots(
        1,
        1,
        clear=True,
    )
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

    ax_m = ax.twinx()
    ax_m.errorbar(
        ji_mb["Time"],
        ji_mb["rate [Gt/yr]"],
        yerr=ji_mb["error [Gt/yr]"],
        color="0",
        linewidth=0.50,
        capsize=1,
        capthick=0.50,
        label="mass loss",
    )
    ax.set_xlabel("Year")
    ax.set_ylabel("Temperature (Celsius)")
    ax.set_xlim(1980, 2021)
    ax.set_ylim(0, 5)
    ax_m.set_ylabel("Dynamic mass loss rate (Gt/yr)")
    l1 = ax.legend(loc="upper left")
    l1.get_frame().set_linewidth(0.0)
    l1.get_frame().set_alpha(0.0)
    l2 = ax_m.legend(loc="upper center")
    l2.get_frame().set_linewidth(0.0)
    l2.get_frame().set_alpha(0.0)
    set_size(3.2, 1.8)

    fig.savefig("jib_obs_1980_2020.pdf", pad_inches="tight")
    ax.hlines(T_1980_1996_mean, 1980, 1997, color="k", linestyle="dashed")
    ax.hlines(T_1997_2015_mean, 1997, 2015, color="k", linestyle="dashed")
    fig.savefig("jib_obs_1980_2020_with_mean.pdf", pad_inches="tight")
