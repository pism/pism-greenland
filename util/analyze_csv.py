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
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
ensemble_file = options.ensemble_file
m_var = options.variable
ifile = options.FILE
calc_second_order = False

df = pd.read_csv(ifile[0], parse_dates=["time"])

if ensemble_file is not None:
    id_df = pd.read_csv(ensemble_file)
    param_names = id_df.drop(columns="id").columns.values.tolist()
    # Define a salib "problem"
    problem = {
        "num_vars": len(id_df.drop(columns="id").columns.values),
        "names": param_names,  # Parameter names
        "bounds": zip(
            id_df.drop(columns="id").min().values, id_df.drop(columns="id").max().values
        ),  # Parameter bounds
    }

    missing_ids = list(set(id_df["id"]).difference(df["id"]))
    if missing_ids:
        print("The following simulation ids are missing:\n   {}".format(missing_ids))

        id_df_missing_removed = id_df[~id_df["id"].isin(missing_ids)]
        id_df_missing = id_df[id_df["id"].isin(missing_ids)]
        params = np.array(
            id_df_missing_removed.drop(columns="id").values, dtype=np.float32
        )
    else:
        params = np.array(id_df.drop(columns="id").values, dtype=np.float32)
        id_df_missing = None

    df = pd.merge(id_df, df, on="id")

m_df = df[
    (df["time"] > pd.to_datetime("1985-1-1"))
    & (df["time"] < pd.to_datetime("1986-1-1"))
]
m_df = df.groupby(by="id").mean().reset_index()
all_ids = m_df["id"].unique()
ids_pass = m_df[m_df[m_var] > -30]["id"]
ids_fail = m_df[~m_df["id"].isin(ids_pass)]

df["pass"] = False
df[df["id"].isin(ids_pass)]["pass"] = True

outside_df = m_df[m_df[m_var] < -30]
inside_df = m_df[m_df[m_var] >= -30]


D = pd.read_csv(
    "~/Google Drive/My Drive/data/mankoff_discharge/gate_merged.csv",
    parse_dates=["Date"],
)
D = D[D["Gate"] == 184]

fig = plt.figure()
ax = fig.add_subplot(111)


def plot_ts(f, ax, m_var="total_grounding_line_flux (Gt year-1)"):
    if f["pass"] is True:

        ax.plot(
            f.index,
            f[m_var],
            color="b",
            lw=0.25,
            alpha=0.2,
        )
    else:
        ax.plot(
            f.index,
            f[m_var],
            color="0.25",
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
ax.set_ylabel("Flux (Gt/yr)")
ax.set_ylim(-75, -10)
set_size(6, 3)
ofile = "calving_{}.pdf".format(m_var)
print("  saving to {}".format(ofile))
fig.savefig(ofile, bbox_inches="tight")


ST_df = []
for s_df in df.groupby(by="time"):
    if id_df_missing is not None:
        response = s_df[1][["id", "total_grounding_line_flux (Gt year-1)"]]
        f = NearestNDInterpolator(params, response.values[:, 1], rescale=True)
        data = f(*np.transpose(id_df_missing.values[:, 1:8]))
        filled = pd.DataFrame(
            data=np.transpose([missing_ids, data]), columns=["id", m_var]
        )
        response_filled = response.append(filled)
        response_filled = response_filled.sort_values(by="id")
        response_matrix = response_filled[response_filled.columns[-1]].values
    else:
        response_matrix = s_df[1]["total_grounding_line_flux (Gt year-1)"].values
    Si = sobol.analyze(
        problem,
        response_matrix,
        calc_second_order=calc_second_order,
        num_resamples=100,
        print_to_console=False,
    )
    if calc_second_order:
        total_Si, first_Si, second_Si = Si.to_df()
    else:
        total_Si, first_Si = Si.to_df()

    t_df = pd.DataFrame(
        data=total_Si["ST"].values.reshape(1, -1), columns=total_Si.transpose().columns
    )
    t_df["Date"] = s_df[0]
    t_df.set_index("Date")
    ST_df.append(t_df)
ST_df = pd.concat(ST_df)
ST_df.reset_index(inplace=True, drop=True)
ST_df.set_index(ST_df["Date"], inplace=True)

fig = plt.figure()
ax = fig.add_subplot(111)
sns.lineplot(data=ST_df, ax=ax).set_title("Sobol Indices")
fig.savefig("total_sobol_indices.pdf", bbox_inches="tight")

fig, axs = plt.subplots(len(param_names), 1, figsize=[4, 12])
fig.subplots_adjust(hspace=0.55, wspace=0.25)
for k, p_var in enumerate(param_names):
    sns.histplot(data=inside_df, x=p_var, stat="count", linewidth=0.8, ax=axs[k])
    ax.set_title(p_var)
    fig.savefig(f"hist_inside_1985.pdf", bbox_inches="tight")

fig, axs = plt.subplots(len(param_names), 1, figsize=[4, 12])
fig.subplots_adjust(hspace=0.55, wspace=0.25)
for k, p_var in enumerate(param_names):
    sns.histplot(data=outside_df, x=p_var, stat="count", linewidth=0.8, ax=axs[k])
    ax.set_title(p_var)
    fig.savefig(f"hist_outside_1985.pdf", bbox_inches="tight")
