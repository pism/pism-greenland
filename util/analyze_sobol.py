#!/usr/bin/env python

# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
import pandas as pd
from pandas.api.types import is_string_dtype
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
parser.add_argument("--second_order", action="store_true", default=False)
parser.add_argument("--variable", default="total_grounding_line_flux (Gt year-1)")
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
calc_second_order = options.second_order
ensemble_file = options.ensemble_file
m_var = options.variable
ifile = options.FILE[0]
m_id = "id"


def prepare_df(ifile):
    df = pd.read_csv(ifile, parse_dates=["time"])

    return df


def compute_sobol_indices(
    df,
    ensemble_file=None,
    calc_second_order=False,
    calc_variable="total_grounding_line_flux (Gt year-1)",
):

    id_df = pd.read_csv(ensemble_file)
    param_names = id_df.drop(columns="id").columns.values.tolist()
    for k, col in id_df.iteritems():
        if is_string_dtype(col):
            u = col.unique()
            u.sort()
            v = [k for k, v in enumerate(u)]
            col.replace(to_replace=dict(zip(u, v)), inplace=True)
    # Define a salib "problem"
    problem = {
        "num_vars": len(id_df.drop(columns="id").columns.values),
        "names": param_names,  # Parameter names
        "bounds": zip(
            id_df.drop(columns="id").min().values,
            id_df.drop(columns="id").max().values,
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

    ST_df = []
    for m_date, s_df in df.groupby(by="time"):
        print(f"Processing {m_date}", s_df.values.shape)
        if id_df_missing is not None:
            response = s_df[["id", m_var]]
            X = id_df_missing.drop(columns="id")
            f = LinearNDInterpolator(params, response.values[:, 1], rescale=True)
            data = f(*np.transpose(X.values))
            filled = pd.DataFrame(
                data=np.transpose([missing_ids, data]), columns=["id", m_var]
            )
            response_filled = response.append(filled)
            response_filled = response_filled.sort_values(by="id")
            response_matrix = response_filled[response_filled.columns[-1]].values
        else:
            response_matrix = s_df[m_var].values
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
            data=total_Si["ST"].values.reshape(1, -1),
            columns=total_Si.transpose().columns,
        )
        t_df["Date"] = m_date
        t_df.set_index("Date")
        ST_df.append(t_df)
    ST_df = pd.concat(ST_df)
    ST_df.reset_index(inplace=True, drop=True)
    ST_df.set_index(ST_df["Date"], inplace=True)
    return ST_df


df = prepare_df(ifile)
ST_df = compute_sobol_indices(
    df,
    ensemble_file=ensemble_file,
    calc_variable=m_var,
    calc_second_order=calc_second_order,
)

fig = plt.figure()
ax = fig.add_subplot(111)
sns.lineplot(data=ST_df, ax=ax).set_title("Sobol Indices")
fig.savefig("total_sobol_indices.pdf", bbox_inches="tight")
