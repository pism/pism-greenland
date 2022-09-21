#!/usr/bin/env python

# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
import pandas as pd
from pandas.api.types import is_string_dtype
import pylab as plt
import seaborn as sns
from SALib.analyze import sobol, delta
import numpy as np
from scipy.interpolate import griddata
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


def prepare_df(ifile):
    df = pd.read_csv(ifile, parse_dates=["time"])

    return df


def compute_sobol_indices(
    df,
    ensemble_file=None,
    calc_second_order=False,
    calc_variable="total_grounding_line_flux (Gt year-1)",
    method="delta",
):

    # remove True/False
    id_df = (pd.read_csv(ensemble_file) * 1).replace(np.nan, 0)

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

    df = pd.merge(id_df, df, on="id")
    Sobol_dfs = []
    for m_date, s_df in df.groupby(by="time"):
        print(f"Processing {m_date}")
        missing_ids = list(set(id_df["id"]).difference(df["id"]))
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
        if id_df_missing is not None:
            response = s_df[["id", m_var]]
            X = id_df_missing.drop(columns="id").values
            data = griddata(params, response.values[:, 1], X, method=interp_method)

            filled = pd.DataFrame(
                data=np.transpose([missing_ids, data]), columns=["id", m_var]
            )
            response_filled = pd.concat([response, filled])
            response_filled = response_filled.sort_values(by="id")
            response_matrix = response_filled[response_filled.columns[-1]].values
        else:
            response_matrix = s_df[m_var].values
        S2_vars = None
        if method == "sobol":
            Si = sobol.analyze(
                problem,
                response_matrix[0:96],
                calc_second_order=calc_second_order,
                num_resamples=100,
                print_to_console=False,
            )
            sobol_indices = ["ST", "S1"]
            Si_df = Si.to_df()
            if calc_second_order:
                sobol_indices.append("S2")
                S2_vars = Si_df[2].index

            s_dfs = []
            for k, s_index in enumerate(sobol_indices):
                m_df = pd.DataFrame(
                    data=Si_df[k][s_index].values.reshape(1, -1),
                    columns=Si_df[k].transpose().columns,
                )
                m_df["Date"] = m_date
                m_df.set_index("Date")
                m_df["Si"] = s_index

                m_conf_df = pd.DataFrame(
                    data=Si_df[k][s_index + "_conf"].values.reshape(1, -1),
                    columns=Si_df[k].transpose().columns,
                )
                m_conf_df["Date"] = m_date
                m_conf_df.set_index("Date")
                m_conf_df["Si"] = s_index + "_conf"
                s_dfs.append(pd.concat([m_df, m_conf_df]))
        elif method == "delta":
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

                m_conf_df = pd.DataFrame(
                    data=Si_df[s_index + "_conf"].values.reshape(1, -1),
                    columns=Si_df.transpose().columns,
                )
                m_conf_df["Date"] = m_date
                m_conf_df.set_index("Date")
                m_conf_df["Si"] = s_index + "_conf"
                s_dfs.append(pd.concat([m_df, m_conf_df]))

        else:
            print(f"Method {method} not implemented")

        s_df = pd.concat(s_dfs)
        Sobol_dfs.append(s_df)

    Sobol_df = pd.concat(Sobol_dfs)
    Sobol_df.reset_index(inplace=True, drop=True)
    Sobol_df.set_index(Sobol_df["Date"], inplace=True)
    return Sobol_df, sobol_indices, S2_vars


# Set up the option parser
parser = ArgumentParser()
parser.description = "A"
parser.add_argument("--ensemble_file", default=None)
parser.add_argument(
    "--interpolation_method", default="nearest", choices=["nearest", "linear"]
)
parser.add_argument("--second_order", action="store_true", default=False)
parser.add_argument("--variable", default="total_grounding_line_flux (Gt year-1)")
parser.add_argument("--method", choices=["sobol", "delta"], default="delta")
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
calc_second_order = options.second_order
ensemble_file = options.ensemble_file
interp_method = options.interpolation_method
method = options.method
m_var = options.variable
ifile = options.FILE[0]
m_id = "id"


df = prepare_df(ifile)
Sobol_df, sobol_indices, S2_vars = compute_sobol_indices(
    df,
    ensemble_file=ensemble_file,
    calc_variable=m_var,
    calc_second_order=calc_second_order,
    method=method,
)
if not calc_second_order:
    S2_vars = ""

fig, axs = plt.subplots(
    len(sobol_indices),
    1,
    sharex="col",
    figsize=[12, 10],
)

for k, si in enumerate(sobol_indices):
    ax = axs[k]
    p_df = (
        Sobol_df[Sobol_df["Si"] == si]
        .drop(columns=S2_vars, errors="ignore")
        .drop(columns=["Date", "Si"])
    )

    p_conf_df = (
        Sobol_df[Sobol_df["Si"] == si + "_conf"]
        .drop(columns=S2_vars, errors="ignore")
        .drop(columns=["Date", "Si"])
    )

    [ax.errorbar(p_df.index, p_df[v], yerr=p_conf_df[v], label=v) for v in p_df.keys()]
if "S2" in sobol_indices:
    ax = axs[-1]
    p_df = Sobol_df[Sobol_df["Si"] == "S2"].drop(columns=["Date", "Si"])
    p_conf_df = Sobol_df[Sobol_df["Si"] == "S2_conf"].drop(columns=["Date", "Si"])
    [ax.errorbar(p_df.index, p_df[v], yerr=p_conf_df[v], label=v) for v in S2_vars]

axs[0].set_title(f"Sobol indices for {m_var}")
legend = axs[0].legend(loc="upper right")
fig.savefig(f"{method}_indices.pdf")
