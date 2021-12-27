#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import uniform, randint, truncnorm, gamma

from SALib.sample import saltelli

parser = ArgumentParser()
parser.description = "Draw samples using the Saltelli methods"
parser.add_argument(
    "-s",
    "--n_samples",
    dest="n_samples",
    type=int,
    help="""number of samples to draw. default=10.""",
    default=10,
)
parser.add_argument(
    "--calc_second_order",
    action="store_true",
    help="""Second order interactions.""",
    default=False,
)
parser.add_argument(
    "OUTFILE",
    nargs=1,
    help="Ouput file (CSV)",
    default="velocity_calibration_samples.csv",
)
options = parser.parse_args()
n_samples = options.n_samples
calc_second_order = options.calc_second_order
outfile = options.OUTFILE[-1]


distributions = {
    "vcm": uniform(loc=0.5, scale=0.75),
    "fracture_softening": uniform(loc=0.40, scale=0.60),
    "fracture_threshold": uniform(loc=40e3, scale=110e3),
    "fracture_healing_rate": uniform(loc=0.0, scale=2.0),
    "fracture_healing_threshold": uniform(loc=1e-11, scale=9.9e-10),
    #  "calving_rate_scaling_file": randint(0, 2),
    #  "frontal_melt_file": randint(0, 10),
    #  "thickness_calving_threshold": randint(200, 400),
    #  "gamma_T": uniform(loc=1.00e-4, scale=0.5e-4),
}


# Generate the Sobol sequence samples with uniform distributions

problem = {
    "num_vars": len(distributions.keys()),
    "names": distributions.keys(),
    "bounds": [[0, 1]] * len(distributions.keys()),
}

# Generate samples
unif_sample = saltelli.sample(problem, n_samples, calc_second_order=calc_second_order)


# To hold the transformed variables
dist_sample = np.zeros_like(unif_sample, dtype="object")


# For each variable, transform with the inverse of the CDF (inv(CDF)=ppf)
for i, key in enumerate(distributions.keys()):
    if key == "calving_rate_scaling_file":
        dist_sample[:, i] = [
            f"seasonal_calving_{int(id)}_1980_2020.nc"
            for id in distributions[key].ppf(unif_sample[:, i])
        ]
    if key == "frontal_melt_file":
        dist_sample[:, i] = [
            f"jib_ocean_forcing_{int(id)}_1980_2020.nc"
            for id in distributions[key].ppf(unif_sample[:, i])
        ]
    else:
        dist_sample[:, i] = distributions[key].ppf(unif_sample[:, i])


# Convert to Pandas dataframe, append column headers, output as csv
df = pd.DataFrame(dist_sample, columns=distributions.keys())
df.to_csv(outfile, index=True, index_label="id")

df["climate"] = "given"
df["hydrology"] = "routing"
df["frontal_melt"] = "discharge_routing"
df["frontal_melt_file"] = "jib_ocean_forcing_ctrl_1980_2020.nc"
df["climate_file"] = "DMI-HIRHAM5_GL2_ERAI_1980_2016_EPSG3413_4500M_DM.nc"
df["runoff_file"] = "DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_4500M_DM.nc"
df["salinity"] = ""
df["pseudo_plastic_q"] = 0.6
df["sia_e"] = 1.25
df["ssa_n"] = 3.0
df["fractures"] = "true"
df["fracture_rate"] = 0.5
df["gamma_T"] = 1.25e-4
df["thickness_calving_threshold"] = 300

df.to_csv(f"ensemble_{outfile}", index=True, index_label="id")
