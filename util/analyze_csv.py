#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype
import pylab as plt
import seaborn as sns
from datetime import datetime


def set_size(w, h, ax=None):
    """w, h: width, height in inches"""

    if not ax:
        ax = plt.gca()
    left = ax.figure.subplotpars.left
    right = ax.figure.subplotpars.right
    top = ax.figure.subplotpars.top
    bottom = ax.figure.subplotpars.bottom
    figw = float(w) / (right - left)
    figh = float(h) / (top - bottom)
    ax.figure.set_size_inches(figw, figh)


def to_decimal_year(date):

    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction


def load_mankoff_discharge():
    """
    Load discharge time-series from Ken Mankoff for
    Jakobshavn Isbrae (Gate 184)
    """

    df = pd.read_csv(
        "~/Google Drive/My Drive/data/mankoff_discharge/gate_merged.csv",
        parse_dates=["Date"],
    )
    df = df[df["Gate"] == 184]
    df["time"] = df["Date"]
    df.set_index("Date", inplace=True)

    return df


def load_mouginot_discharge():
    """
    Load discharge time-series from Jeremy Mouginot
    """
    mou19_d = pd.read_excel(
        "/Users/andy/Google Drive/My Drive/data/mouginot_discharge/pnas.1904242116.sd02.xlsx",
        sheet_name="(5) year_D_R23p2-5.5km",
        header=2,
        usecols="A,P:BJ",
        engine="openpyxl",
    )
    mou19_d_err = pd.read_excel(
        "/Users/andy/Google Drive/My Drive/data/mouginot_discharge/pnas.1904242116.sd02.xlsx",
        sheet_name="(5) year_D_R23p2-5.5km",
        header=2,
        usecols="A,BN:DH",
        engine="openpyxl",
    )

    mou19_d = mou19_d[mou19_d["NAME"] == "JAKOBSHAVN_ISBRAE"].drop(columns={"NAME"})
    mou19_d_err = mou19_d_err[mou19_d_err["NAME"] == "JAKOBSHAVN_ISBRAE"].drop(
        columns={"NAME"}
    )
    years = mou19_d.columns.astype("float").values - 0.5
    d = mou19_d.values.T
    d_err = mou19_d_err.values.T

    mou19_d = pd.DataFrame(
        data=d, index=years.astype(int), columns={"Discharge (Gt/yr)"}
    )
    mou19_d_err = pd.DataFrame(
        data=d_err, index=years.astype(int), columns={"Discharge Error (Gt/yr)"}
    )

    df = pd.merge(mou19_d, mou19_d_err, left_index=True, right_index=True)
    dates = pd.to_datetime(df.index, format="%Y")
    df.set_index(dates, inplace=True)
    df["time"] = dates

    return df


def load_data(data_file, ensemble_file=None):

    df = pd.read_csv(data_file, parse_dates=["time"])
    df[m_var] *= -1

    id_df_missing = None
    if ensemble_file is not None:
        id_df = pd.read_csv(ensemble_file)
        param_names = id_df.drop(columns="id").columns.values.tolist()
        for k, col in id_df.iteritems():
            if is_string_dtype(col):
                u = col.unique()
                u.sort()
                v = [k for k, v in enumerate(u)]
                col.replace(to_replace=dict(zip(u, v)), inplace=True)
        # Define a salib "problem"

        missing_ids = list(set(id_df["id"]).difference(df["id"]))
        if missing_ids:
            print(
                "The following simulation ids are missing:\n   {}".format(missing_ids)
            )

            id_df_missing_removed = id_df[~id_df["id"].isin(missing_ids)]
            id_df_missing = id_df[id_df["id"].isin(missing_ids)]

        df = pd.merge(id_df, df, on="id")
    return df, param_names


if __name__ == "__main__":

    # Set up the option parser
    parser = ArgumentParser()
    parser.description = "A"
    parser.add_argument("--beta", default=1.0)
    parser.add_argument(
        "--ensemble_file",
        default="../historical/2022_01_fractures_melt/uq/jib_fractures_melt.csv",
    )
    parser.add_argument("--variable", default="total_grounding_line_flux (Gt year-1)")
    parser.add_argument(
        "--mean_file",
        default="../historical/2022_01_fractures_melt/csv/fldmean_ts.csv",
    )
    parser.add_argument(
        "--sum_file", default="../historical/2022_01_fractures_melt/csv/fldsum_ts.csv"
    )
    parser.add_argument("--threshold_range", default=[20, 40])
    parser.add_argument("--smoothing_length", type=int, default=1)
    options = parser.parse_args()
    beta = options.beta
    ensemble_file = options.ensemble_file
    m_var = options.variable
    calc_second_order = False
    smoothing_length = options.smoothing_length
    threshold = options.threshold_range

    d_man = load_mankoff_discharge()
    d_mou = load_mouginot_discharge()

    mean_df, _ = load_data(options.mean_file, ensemble_file=ensemble_file)
    sum_df, param_names = load_data(options.sum_file, ensemble_file=ensemble_file)

    validation_time = [pd.to_datetime("1985-1-1"), pd.to_datetime("1986-1-1")]
    m_df = sum_df[sum_df["time"].between(*validation_time)]
    m_df = sum_df.groupby(by="id").mean().reset_index()
    all_ids = sum_df["id"].unique()
    d_val = d_mou[d_mou["time"].between(*validation_time)]
    discharge_range = np.array(
        [
            (d_val["Discharge (Gt/yr)"].max() - d_val["Discharge Error (Gt/yr)"].max()),
            (d_val["Discharge (Gt/yr)"].max() + d_val["Discharge Error (Gt/yr)"].max()),
        ]
    )
    discharge_range *= beta
    ids_pass = m_df[m_df[m_var].between(*discharge_range)]["id"]
    ids_fail = m_df[~m_df["id"].isin(ids_pass)]["id"]

    sum_df["pass"] = False
    sum_df["pass"].loc[sum_df["id"].isin(ids_pass)] = True

    fig = plt.figure(
        figsize=[6.0, 2.2],
    )
    ax = fig.add_subplot(111)

    def plot_ts(df, ax, m_var, **kwargs):
        ax.plot(
            df.index,
            df[m_var],
            smoothing_length,
            alpha=0.5,
            **kwargs,
        )

    if smoothing_length > 1:
        kwargs = {"color": "0.25", "lw": 0.25}
        [
            plot_ts(f, ax, m_var, **kwargs)
            for f in sum_df.groupby(by="id").rolling(smoothing_length, on="time")
        ]
        kwargs = {"color": "#6a51a3", "lw": 0.5}
        [
            plot_ts(f, ax, m_var, **kwargs)
            for f in sum_df[mean_df["id"].isin(ids_pass)]
            .groupby(by="id")
            .rolling(smoothing_length, on="time")
        ]
    else:
        kwargs = {"color": ".25", "lw": 0.25}
        [
            plot_ts(f[1].set_index("time"), ax, m_var, **kwargs)
            for f in sum_df.groupby(by="id")
        ]
        kwargs = {"color": "#6a51a3", "lw": 0.5}
        [
            plot_ts(f[1].set_index("time"), ax, m_var, **kwargs)
            for f in sum_df[sum_df["id"].isin(ids_pass)].groupby(by="id")
        ]

    ax.fill_between(
        d_man.index,
        d_man["Discharge [Gt/yr]"] - d_man["Discharge Error [Gt/yr]"],
        d_man["Discharge [Gt/yr]"] + d_man["Discharge Error [Gt/yr]"],
        color="#238b45",
        alpha=0.75,
        lw=1.0,
    )

    ax.fill_between(
        d_mou.index,
        d_mou["Discharge (Gt/yr)"] - d_mou["Discharge Error (Gt/yr)"],
        d_mou["Discharge (Gt/yr)"] + d_mou["Discharge Error (Gt/yr)"],
        color="#bdd7e7",
        alpha=0.75,
        lw=1.0,
    )

    ax.plot(
        d_man.index,
        d_man["Discharge [Gt/yr]"],
        color="#238b45",
        lw=2.0,
        label="D (Mankoff)",
    )
    ax.plot(
        d_mou.index,
        d_mou["Discharge (Gt/yr)"],
        color="#2171b5",
        lw=2.0,
        label="D (Mouginot)",
    )

    ax.set_xlim(datetime(1980, 1, 1), datetime(2020, 1, 1))
    ax.set_xlabel("Year")
    ax.set_ylabel("Flux (Gt/yr)")
    ax.set_ylim(0, 100)
    legend = ax.legend()
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.0)

    set_size(6, 3)
    ofile = "jib_{}.pdf".format(m_var.split(" ")[0])
    print("  saving to {}".format(ofile))
    fig.savefig(ofile, bbox_inches="tight")

    fig = plt.figure(
        figsize=[6.0, 2.2],
    )
    ax = fig.add_subplot(111)

    m_var = "vonmises_calving_rate (m year-1)"

    if smoothing_length > 1:
        kwargs = {"color": "0.25", "lw": 0.25}
        [
            plot_ts(f, ax, m_var, **kwargs)
            for f in mean_df.groupby(by="id").rolling(smoothing_length, on="time")
        ]
        kwargs = {"color": "#6a51a3", "lw": 0.5}
        [
            plot_ts(f, ax, m_var, color=color)
            for f in mean_df[mean_df["id"].isin(ids_pass)]
            .groupby(by="id")
            .rolling(smoothing_length, on="time")
        ]
    else:
        kwargs = {"color": "0.25", "lw": 0.25}
        [
            plot_ts(f[1].set_index("time"), ax, m_var, **kwargs)
            for f in mean_df.groupby(by="id")
        ]
        kwargs = {"color": "#6a51a3", "lw": 0.5}
        [
            plot_ts(f[1].set_index("time"), ax, m_var, **kwargs)
            for f in mean_df[mean_df["id"].isin(ids_pass)].groupby(by="id")
        ]

    ax.set_xlim(datetime(1980, 1, 1), datetime(1990, 1, 1))
    ax.set_xlabel("Year")
    ax.set_ylabel("Calving rate (m/yr)")
    # ax.set_ylim(0, 100)
    set_size(6, 3)
    ofile = "jib_{}.pdf".format(m_var.split(" ")[0])
    print("  saving to {}".format(ofile))
    fig.savefig(ofile, bbox_inches="tight")

    fig, axs = plt.subplots(len(param_names), 1, figsize=[4, 15])
    fig.subplots_adjust(hspace=0.55, wspace=0.25)
    for k, p_var in enumerate(param_names):
        sns.histplot(
            data=sum_df,
            x=p_var,
            hue="pass",
            multiple="layer",
            linewidth=0.8,
            ax=axs[k],
            legend=False,
        )
        ax.set_title(p_var)
        fig.savefig(f"params_hist_1985.pdf", bbox_inches="tight")
