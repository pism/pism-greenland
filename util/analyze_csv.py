#!/usr/bin/env python

# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser
import pandas as pd
import pylab as plt
import seaborn as sns
from SALib.analyze import sobol
import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator


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
    "names": id_df.columns[1:8],  # Parameter names
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


for s_df in df.groupby(by="time"):
    response = s_df[1][["id", "total_grounding_line_flux (Gt year-1)"]]
    f = NearestNDInterpolator(params, response.values[:, 1], rescale=True)
    data = f(*np.transpose(id_df_missing.values[:, 1:8]))
    filled = pd.DataFrame(data=np.transpose([missing_ids, data]), columns=["id", m_var])
    response_filled = response.append(filled)
    response_filled = response_filled.sort_values(by="id")
    response_matrix = response_filled[response_filled.columns[-1]].values

    Si = sobol.analyze(problem, response_matrix, calc_second_order=True, num_resamples=100, print_to_console=False)

# for p_var in [
#     "vcm",
#     "fracture_softening",
#     "fracture_rate",
#     "fracture_threshold",
#     "fracture_healing_rate",
#     "fracture_healing_threshold",
# ]:
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     sns.histplot(data=outside_df, p=m_var, stat="density", linewidth=0.8, ax=ax)
