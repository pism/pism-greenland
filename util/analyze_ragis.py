#!/bin/env python3

# Copyright (C) 2023 Andy Aschwanden
#
# This file is part of pism-greenland.
#
# PISM-EMULATOR is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM-EMULATOR is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import numpy as np
import pylab as plt
from glob import glob
import pandas as pd
from pandas.api.types import is_string_dtype
import os
import re
from scipy.interpolate import interp1d
import seaborn as sns
from SALib.analyze import delta
from sklearn.metrics import mean_squared_error
import xarray as xr
import xcdat
from datetime import datetime
from argparse import ArgumentParser


def to_decimal_year(date):
    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction


def load_imbie():
    """
    Loading the IMBIE Greenland data set downloaded from
    http://imbie.org/wp-content/uploads/2012/11/imbie_dataset_greenland_dynamics-2020_02_28.xlsx

    """
    df_df = pd.read_excel(
        "http://imbie.org/wp-content/uploads/2012/11/imbie_dataset_greenland_dynamics-2020_02_28.xlsx",
        sheet_name="Greenland Ice Mass",
        engine="openpyxl",
    )
    df = df_df[
        [
            "Year",
            "Cumulative ice sheet mass change (Gt)",
            "Cumulative ice sheet mass change uncertainty (Gt)",
            "Cumulative surface mass balance anomaly (Gt)",
            "Cumulative surface mass balance anomaly uncertainty (Gt)",
            "Cumulative ice dynamics anomaly (Gt)",
            "Cumulative ice dynamics anomaly uncertainty (Gt)",
            "Rate of mass balance anomaly (Gt/yr)",
            "Rate of ice dynamics anomaly (Gt/yr)",
            "Rate of mass balance anomaly uncertainty (Gt/yr)",
            "Rate of ice dyanamics anomaly uncertainty (Gt/yr)",
        ]
    ].rename(
        columns={
            "Cumulative ice sheet mass change (Gt)": "Mass (Gt)",
            "Cumulative ice sheet mass change uncertainty (Gt)": "Mass uncertainty (Gt)",
            "Rate of mass balance anomaly (Gt/yr)": "SMB (Gt/yr)",
            "Rate of ice dynamics anomaly (Gt/yr)": "D (Gt/yr)",
            "Rate of mass balance anomaly uncertainty (Gt/yr)": "SMB uncertainty (Gt/yr)",
            "Rate of ice dyanamics anomaly uncertainty (Gt/yr)": "D uncertainty (Gt/yr)",
        }
    )

    df = df[df["Year"] >= 1992.0]
    s = df[(df["Year"] >= 1980) & (df["Year"] < 1990)]
    mass_mean = s["Mass (Gt)"].mean() / (1990 - 1980)
    smb_mean = s["Cumulative surface mass balance anomaly (Gt)"].mean() / (1990 - 1980)
    df[f"SMB (Gt/yr)"] += 2 * 1964 / 10
    df[f"D (Gt/yr)"] -= 2 * 1964 / 10
    cmSLE = 1.0 / 362.5 / 10.0
    df[f"SLE (cm)"] = -df["Mass (Gt)"] * cmSLE
    df[f"SLE uncertainty (cm)"] = df["Mass uncertainty (Gt)"] * cmSLE

    y = df["Year"].astype("int")
    df["Date"] = pd.to_datetime({"year": y, "month": 1, "day": 1}) + pd.to_timedelta(
        (df["Year"] - df["Year"].astype("int")) * 3.15569259747e7, "seconds"
    )

    return df


def normalize_by_date(ds, varname="limnsw", date="1992-1-1"):
    return (
        ds[varname] - ds.sel(time=pd.to_datetime("1992-1-1"), method="nearest")[varname]
    )


if __name__ == "__main__":
    __spec__ = None

    parser = ArgumentParser()
    parser.add_argument("FILES", nargs="*", help="Time series files", default=None)
    parser.add_argument(
        "-o", dest="outfile", help="Output filename", default="ragis_scalar_ts.pdf"
    )

    args = parser.parse_args()
    ts_files = sorted(args.FILES)
    outfile = args.outfile

    mass_varname = "limnsw"
    flux_varname = "grounding_line_flux"
    bg_color = "#216779"
    imbie_color = "#84cedb"

    imbie = load_imbie()

    kg2cmsle = 1 / 1e12 * 1.0 / 362.5 / 10.0
    gt2cmsle = 1 / 362.5 / 10.0
    sigma = 2
    plt.style.use("dark_background")

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex="col", figsize=(12, 8))
    fig.subplots_adjust(wspace=0.0, hspace=0.0)

    labels = []
    dfs = []
    for m_file in ts_files:
        print(f"Reading {m_file}")
        with xcdat.open_dataset(m_file) as ds:
            mass_var = normalize_by_date(ds) / 1e12
            mass_df = mass_var.to_dataframe().rename(columns={"limnsw": "Mass (Gt)"})
            sle_var = mass_var * gt2cmsle
            d_df = (
                ds[flux_varname]
                .to_dataframe()
                .rename(columns={"grounding_line_flux": "D (Gt/yr)"})
            )
            df = pd.merge(mass_df, d_df, on="time")
            years = [to_decimal_year(x) for x in mass_df.index]
            df["Year"] = years

            if hasattr(ds, "id"):
                m_id = ds.attrs["id"]
            else:
                m_id = int(re.search("id_(.+?)_", m_file).group(1))
            try:
                m_id = int(m_id)
            except:
                pass
            df["Experiment"] = m_id
            dfs.append(df)
            l = sle_var.plot(ax=axs[0], color="w", lw=0.5)[0]

            labels.append(l.get_label())
            flux_var_running = (
                ds[flux_varname].rolling(time=390, center=True).mean().dropna("time")
            )
            flux_var_running.plot(ax=axs[1], color="w", lw=0.5)

    axs[0].fill_between(
        imbie["Date"],
        (imbie["Mass (Gt)"] + sigma * imbie["Mass uncertainty (Gt)"]) * gt2cmsle,
        (imbie["Mass (Gt)"] - sigma * imbie["Mass uncertainty (Gt)"]) * gt2cmsle,
        ls="solid",
        lw=0,
        alpha=1,
        color=imbie_color,
        label="2-$\sigma$ IMBIE",
    )
    axs[1].fill_between(
        imbie["Date"],
        (imbie["D (Gt/yr)"] + sigma * imbie["D uncertainty (Gt/yr)"]),
        (imbie["D (Gt/yr)"] - sigma * imbie["D uncertainty (Gt/yr)"]),
        ls="solid",
        lw=0,
        alpha=1,
        color=imbie_color,
    )

    fig.set_facecolor(bg_color)
    for ax in axs:
        ax.set_facecolor(bg_color)
    fig.set_facecolor(bg_color)

    axs[0].set_xlabel("")
    axs[0].set_ylabel("Contribution to sea-level since 1992\n(cm SLE)")
    axs[1].set_xlabel("Year")
    axs[1].set_ylabel("Solid Discharge (Gt/yr)")
    axs[0].set_xlim(
        pd.to_datetime("1980-1-1", format="%Y-%m-%d"),
        pd.to_datetime("2020-1-1", format="%Y-%m-%d"),
    )
    axs[0].set_ylim(-1.5, 1.5)
    axs[1].set_ylim(-3000, 200)
    legend = axs[0].legend()
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.0)
    fig.savefig(outfile)

    data_df = pd.concat(dfs)
    lenghts = []
    exps = []
    for exp, df in data_df.groupby(by="Experiment"):
        exps.append(exp)
        lenghts.append(len(df))
    max_length = max(lenghts)
    successful_runs = (np.array(lenghts) == max_length).nonzero()[0]
    data_df = data_df[data_df["Experiment"].isin(successful_runs)]
