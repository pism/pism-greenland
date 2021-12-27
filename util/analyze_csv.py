#!/usr/bin/env python

# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
import pandas as pd
import pylab as plt
import seaborn as sns
from SALib.analyze import sobol
import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
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


def to_decimal_year(date):

    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction


# Set up the option parser
parser = ArgumentParser()
parser.description = "A"
parser.add_argument("--ensemble_file", default=None)
parser.add_argument("--variable", default="total_grounding_line_flux (Gt year-1)")
parser.add_argument("--smoothing_length", type=int, default=1)
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
ensemble_file = options.ensemble_file
m_var = options.variable
ifile = options.FILE
calc_second_order = False
smoothing_length = options.smoothing_length

df = pd.read_csv(ifile[0], parse_dates=["time"])

id_df_missing = None
if ensemble_file is not None:
    id_df = pd.read_csv(ensemble_file)
    param_names = id_df.drop(columns="id").columns.values.tolist()
    # Define a salib "problem"

    missing_ids = list(set(id_df["id"]).difference(df["id"]))
    if missing_ids:
        print("The following simulation ids are missing:\n   {}".format(missing_ids))

        id_df_missing_removed = id_df[~id_df["id"].isin(missing_ids)]
        id_df_missing = id_df[id_df["id"].isin(missing_ids)]

    df = pd.merge(id_df, df, on="id")

m_df = df[
    (df["time"] > pd.to_datetime("1985-1-1"))
    & (df["time"] < pd.to_datetime("1986-1-1"))
]
m_df = df.groupby(by="id").mean().reset_index()
all_ids = m_df["id"].unique()
ids_pass = m_df[m_df[m_var] > -50]["id"]
ids_fail = m_df[~m_df["id"].isin(ids_pass)]

df["pass"] = False
df[df["id"].isin(ids_pass)]["pass"] = True

outside_df = m_df[m_df[m_var] < -40]
inside_df = m_df[m_df[m_var] >= -40]


D = pd.read_csv(
    "~/Google Drive/My Drive/data/mankoff_discharge/gate_merged.csv",
    parse_dates=["Date"],
)
D = D[D["Gate"] == 184]

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

mou19_D = pd.DataFrame(data=d, index=years.astype(int), columns={"Discharge (Gt/yr)"})
mou19_D_err = pd.DataFrame(
    data=d_err, index=years.astype(int), columns={"Discharge Error (Gt/yr)"}
)

mou19 = pd.merge(mou19_D, mou19_D_err, left_index=True, right_index=True)
mou19_dates = pd.to_datetime(mou19.index, format="%Y")
mou19.set_index(mou19_dates, inplace=True)

fig = plt.figure()
ax = fig.add_subplot(111)


def plot_ts(df, ax, m_var):
    ax.plot(
        df.index,
        -df[m_var],
        color="0.25",
        lw=0.5,
        alpha=0.5,
    )


if smoothing_length > 1:
    [
        plot_ts(f[1], ax, m_var)
        for f in df.groupby(by="id").rolling(smoothing_length, on="time")
    ]
else:
    [plot_ts(f[1].set_index("time"), ax, m_var) for f in df.groupby(by="id")]

ax.fill_between(
    D["Date"],
    -D["Discharge [Gt/yr]"] - D["Discharge Error [Gt/yr]"],
    -D["Discharge [Gt/yr]"] + D["Discharge Error [Gt/yr]"],
    color="#238b45",
    alpha=0.5,
    lw=1.0,
)

ax.fill_between(
    mou19.index,
    mou19["Discharge (Gt/yr)"] - mou19["Discharge Error (Gt/yr)"],
    mou19["Discharge (Gt/yr)"] + mou19["Discharge Error (Gt/yr)"],
    color="#bdd7e7",
    alpha=0.5,
    lw=1.0,
)

ax.plot(
    D["Date"],
    D["Discharge [Gt/yr]"],
    color="#238b45",
    lw=1.0,
    label="Mankoff",
)
ax.plot(
    mou19.index,
    mou19["Discharge (Gt/yr)"],
    color="#2171b5",
    lw=1.0,
    label="Mouginot",
)

ax.set_xlim(datetime(1980, 1, 1), datetime(2000, 1, 1))
ax.set_xlabel("Year")
ax.set_ylabel("Flux (Gt/yr)")
ax.set_ylim(0, 75)
set_size(6, 3)
ofile = "calving_{}.pdf".format(m_var)
print("  saving to {}".format(ofile))
fig.savefig(ofile, bbox_inches="tight")


fig, axs = plt.subplots(len(param_names), 1, figsize=[4, 12])
fig.subplots_adjust(hspace=0.55, wspace=0.25)
for k, p_var in enumerate(param_names):
    sns.histplot(data=outside_df, x=p_var, stat="count", linewidth=0.8, ax=axs[k])
    ax.set_title(p_var)
    fig.savefig(f"hist_outside_1985.pdf", bbox_inches="tight")
