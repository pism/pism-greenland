#!/usr/bin/env python

# Copyright (C) 2021-23 Andy Aschwanden

from argparse import ArgumentParser
import pandas as pd
from pandas.api.types import is_string_dtype
import pylab as plt
import seaborn as sns
from SALib.analyze import sobol, delta
import numpy as np
from scipy.interpolate import griddata
from datetime import datetime
import pathlib


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


def prepare_df(ifile):
    suffix = pathlib.Path(ifile).suffix
    if suffix in (".csv", ".gz"):
        df = pd.read_csv(ifile, parse_dates=["time"])
    elif suffix in (".parquet"):
        df = pd.read_parquet(ifile)
    else:
        print(f"{suffix} not recognized")

    return df


def compute_sobol_indices(
    df,
    ensemble_file=None,
    calc_variables=[
        "grounding_line_flux (Gt year-1)",
        "tendency_of_ice_mass_due_to_calving (Gt year-1)",
    ],
):

    # remove True/False
    id_df = (pd.read_csv(ensemble_file) * 1).replace(np.nan, 0)

    param_names = id_df.drop(columns="id").columns.values.tolist()
    for k, col in id_df.items():
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

    df = pd.merge(id_df, df, on="id")
    Sobol_dfs = []
    for m_date, s_df in df.groupby(by="time"):
        print(f"Processing {m_date}")
        missing_ids = list(set(id_df["id"]).difference(s_df["id"]))
        if missing_ids:
            print(
                "The following simulation ids are missing:\n   {}".format(missing_ids)
            )

            id_df_missing_removed = id_df[~id_df["id"].isin(missing_ids)]
            id_df_missing = id_df[id_df["id"].isin(missing_ids)]
            params = np.array(
                id_df_missing_removed.drop(columns="id").values, dtype=np.float32
            )
        else:
            params = np.array(id_df.drop(columns="id").values, dtype=np.float32)
            id_df_missing = None
        for calc_variable in calc_variables:
            if id_df_missing is not None:
                if method == "sobol":
                    response = s_df[["id", calc_variable]]
                    X = id_df_missing.drop(columns="id").values
                    data = griddata(
                        params, response.values[:, 1], X, method=interp_method
                    )

                    filled = pd.DataFrame(
                        data=np.transpose([missing_ids, data]),
                        columns=["id", calc_variable],
                    )
                    response_filled = pd.concat([response, filled])
                    response_filled = response_filled.sort_values(by="id")
                    response_matrix = response_filled[
                        response_filled.columns[-1]
                    ].values
                else:
                    id_df = id_df_missing_removed
                    response_matrix = s_df[calc_variable].values
            else:
                response_matrix = s_df[calc_variable].values
            Si = delta.analyze(
                problem,
                id_df.drop(columns=["id"]).values,
                response_matrix,
                num_resamples=100,
                print_to_console=False,
            )
            sobol_indices = ["delta", "S1"]
            Si_df = Si.to_df()

            s_dfs = []
            for s_index in sobol_indices:
                m_df = pd.DataFrame(
                    data=Si_df[s_index].values.reshape(1, -1),
                    columns=Si_df.transpose().columns,
                )
                m_df["Date"] = m_date
                m_df.set_index("Date")
                m_df["Si"] = s_index
                m_df["Variable"] = calc_variable

                m_conf_df = pd.DataFrame(
                    data=Si_df[s_index + "_conf"].values.reshape(1, -1),
                    columns=Si_df.transpose().columns,
                )
                m_conf_df["Date"] = m_date
                m_conf_df.set_index("Date")
                m_conf_df["Si"] = s_index + "_conf"
                m_conf_df["Variable"] = calc_variable
                s_dfs.append(pd.concat([m_df, m_conf_df]))

            a_df = pd.concat(s_dfs)
        Sobol_dfs.append(a_df)

    Sobol_df = pd.concat(Sobol_dfs)
    Sobol_df.reset_index(inplace=True, drop=True)
    Sobol_df.set_index(Sobol_df["Date"], inplace=True)
    return Sobol_df, sobol_indices


# Set up the option parser
parser = ArgumentParser()
parser.description = "A"
parser.add_argument("--ensemble_file", default=None)
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
ensemble_file = options.ensemble_file
ifile = options.FILE[0]
m_id = "id"


df = prepare_df(ifile)
calc_variables = df.drop(columns=["time", "id"]).columns

Sobol_df, sobol_indices = compute_sobol_indices(
    df,
    ensemble_file=ensemble_file,
)

fig, axs = plt.subplots(
    len(sobol_indices),
    1,
    sharex="col",
    figsize=[12, 10],
)

for k, si in enumerate(sobol_indices):
    m_df = Sobol_df[Sobol_df["Variable"] == "grounding_line_flux (Gt year-1)")]
    ax = axs[k]
    p_df = m_df[m_df["Si"] == si].drop(columns=["Date", "Si"])

    p_conf_df = m_df[m_df["Si"] == si + "_conf"].drop(columns=["Date", "Si"])

    [ax.errorbar(p_df.index, p_df[v], yerr=p_conf_df[v], label=v) for v in p_df.keys()]

axs[0].set_title(f"Sobol indices for {m_var}")
legend = axs[0].legend(loc="upper right")
fig.savefig(f"{method}_indices.pdf")
