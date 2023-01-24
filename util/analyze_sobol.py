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
from joblib import Parallel, delayed

import contextlib
import joblib
from tqdm import tqdm


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


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


def run_analysis(
    df,
    ensemble_file=None,
    calc_variables=["grounding_line_flux (Gt year-1)", "limnsw (kg)"],
    n_jobs=4,
    sobol_indices=["delta", "S1"],
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
    # filter out dates with only 1 experiment, e.g., due to
    # CDO changing the mid-point time when running averaging
    df = pd.concat([x for _, x in df.groupby(by="time") if len(x) > 1])
    n_dates = len(df["time"].unique())
    if n_jobs == 1:
        Sobol_dfs = []
        for m_date, s_df in df.groupby(by="time"):
            Sobol_dfs.append(
                compute_sobol_indices(
                    m_date,
                    s_df,
                    id_df,
                    problem,
                    calc_variables,
                    sobol_indices=sobol_indices,
                )
            )
    else:
        with tqdm_joblib(tqdm(desc="Processing file", total=n_dates)) as progress_bar:
            Sobol_dfs = Parallel(n_jobs=n_jobs)(
                delayed(compute_sobol_indices)(
                    m_date,
                    s_df,
                    id_df,
                    problem,
                    calc_variables,
                    sobol_indices=sobol_indices,
                )
                for m_date, s_df in df.groupby(by="time")
            )

    Sobol_df = pd.concat(Sobol_dfs)
    Sobol_df.reset_index(inplace=True, drop=True)
    return Sobol_df, sobol_indices


def compute_sobol_indices(
    m_date,
    s_df,
    id_df,
    problem,
    calc_variables,
    sobol_indices=["delta", "S1"],
):
    print(f"Processing {m_date}")
    missing_ids = list(set(id_df["id"]).difference(s_df["id"]))
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
    Sobol_dfs = []
    for calc_variable in calc_variables:
        response_matrix = s_df[calc_variable].values
        Si = delta.analyze(
            problem,
            params,
            response_matrix,
            num_resamples=100,
            print_to_console=False,
        )
        Si_df = Si.to_df()

        s_dfs = []
        for s_index in sobol_indices:
            m_df = pd.DataFrame(
                data=Si_df[s_index].values.reshape(1, -1),
                columns=Si_df.transpose().columns,
            )
            m_df["Date"] = m_date
            m_df["Si"] = s_index
            m_df["Variable"] = calc_variable

            m_conf_df = pd.DataFrame(
                data=Si_df[s_index + "_conf"].values.reshape(1, -1),
                columns=Si_df.transpose().columns,
            )
            m_conf_df["Date"] = m_date
            m_conf_df["Si"] = s_index + "_conf"
            m_conf_df["Variable"] = calc_variable
            s_dfs.append(pd.concat([m_df, m_conf_df]))

        a_df = pd.concat(s_dfs)
        Sobol_dfs.append(a_df)
    return pd.concat(Sobol_dfs)


# Set up the option parser
parser = ArgumentParser()
parser.description = "A"
parser.add_argument("--ensemble_file", default=None)
parser.add_argument("-n", "--n_jobs", type=int, default=4)
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
ensemble_file = options.ensemble_file
ifile = options.FILE[0]
n_jobs = options.n_jobs
m_id = "id"


df = prepare_df(ifile)
calc_variables = df.drop(columns=["time", "id"]).columns

Sobol_df, sobol_indices = run_analysis(df, ensemble_file=ensemble_file, n_jobs=n_jobs)
si = "delta"

plt.style.use("fivethirtyeight")


fig, axs = plt.subplots(
    2,
    1,
    sharex="col",
    figsize=[12, 10],
)
fig.subplots_adjust(bottom=0.16)
for k, m_var in enumerate(["limnsw (kg)", "grounding_line_flux (Gt year-1)"]):
    m_df = Sobol_df[Sobol_df["Variable"] == m_var]
    ax = axs.ravel()[k]
    p_df = m_df[m_df["Si"] == si].drop(columns=["Si", "Variable"]).set_index("Date")
    p_conf_df = m_df[m_df["Si"] == si + "_conf"].drop(columns=["Si"])

    [
        ax.errorbar(p_df.index, p_df[v], yerr=p_conf_df[v], lw=2, label=v)
        for v in [
            "vcm",
            "gamma_T",
            "thickness_calving_threshold",
            "ocean_file",
            "sia_e",
            "ssa_n",
            "pseudo_plastic_q",
            "till_effective_fraction_overburden",
        ]
    ]

    lgd = ax.set_title(f"{si} indices for '{m_var}'")
legend = axs[-1].legend(loc="lower left", ncols=3, bbox_to_anchor=(0, -0.4))
fig.tight_layout()
fig.savefig(f"{si}_indices.pdf")
