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
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
ensemble_file = options.ensemble_file
ifile = options.FILE
m_var = "total_grounding_line_flux (Gt year-1)"

df = pd.read_csv(ifile[0], parse_dates=["time"])
id_df = pd.read_csv(ensemble_file)

# Define a salib "problem"
problem = {
    "num_vars": len(id_df.columns[1:8].values),
    "names": id_df.columns[1:8].values.tolist(),  # Parameter names
    "bounds": zip(id_df.min()[1:8], id_df.max()[1:8]),  # Parameter bounds
}


missing_ids = list(set(id_df["id"]).difference(df["id"]))
if missing_ids:
    print("The following simulation ids are missing:\n   {}".format(missing_ids))

    id_df_missing_removed = id_df[~id_df["id"].isin(missing_ids)]
    id_df_missing = id_df[id_df["id"].isin(missing_ids)]
    params = np.array(id_df_missing_removed.values[:, 1:8], dtype=np.float32)


df = pd.merge(id_df, df, on="id")

m_df = df[(df["time"] > pd.to_datetime("1985-1-1")) & (df["time"] < pd.to_datetime("1986-1-1"))]
m_df = df.groupby(by="id").mean().reset_index()

outside_df = m_df[m_df[m_var] < -30]

ST_df = []
for s_df in df.groupby(by="time"):
    response = s_df[1][["id", "total_grounding_line_flux (Gt year-1)"]]
    f = NearestNDInterpolator(params, response.values[:, 1], rescale=True)
    data = f(*np.transpose(id_df_missing.values[:, 1:8]))
    filled = pd.DataFrame(data=np.transpose([missing_ids, data]), columns=["id", m_var])
    response_filled = response.append(filled)
    response_filled = response_filled.sort_values(by="id")
    response_matrix = response_filled[response_filled.columns[-1]].values

    Si = sobol.analyze(problem, response_matrix, calc_second_order=True, num_resamples=100, print_to_console=False)
    total_Si, first_Si, second_Si = Si.to_df()
    t_df = pd.DataFrame(data=total_Si["ST"].values.reshape(1, -1), columns=total_Si.transpose().columns)
    t_df["date"] = s_df[0]
    t_df.set_index("date")
    ST_df.append(t_df)
ST_df = pd.concat(ST_df)
ST_df.reset_index(inplace=True, drop=True)
time = pd.date_range(start="01-15-1980", end="01-01-1986", freq="M")
ST_df["time"] = time
ST_df.set_index(time, inplace=True)

fig = plt.figure()
ax = fig.add_subplot(111)
sns.lineplot(data=ST_df, ax=ax).set_title("Sobol Indices")
fig.savefig("total_sobol_indices.pdf", bbox_inches="tight")

fig, axs = plt.subplots(6, 1, figsize=[4, 12])
fig.subplots_adjust(hspace=0.55, wspace=0.25)
for k, p_var in enumerate(
    [
        "vcm",
        "fracture_softening",
        "fracture_rate",
        "fracture_threshold",
        "fracture_healing_rate",
        "fracture_healing_threshold",
    ]
):
    sns.histplot(data=outside_df, x=p_var, stat="count", linewidth=0.8, ax=axs[k])
    ax.set_title(p_var)
    fig.savefig(f"hist_1985.pdf", bbox_inches="tight")
