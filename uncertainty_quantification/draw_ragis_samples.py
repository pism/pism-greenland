#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import uniform, randint, truncnorm, gamma

from SALib.sample import saltelli


dists = {
    "calving": {
        "uq": {
            "vcm": uniform(loc=0.25, scale=0.75),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "ocean_file": "MAR3.9_MIROC-ESM-CHEM_rcp85_ocean_1960-2100_v4.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_MM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_MM.nc",
            "gamma_T": 1.00e-4,
            "salinity": 34,
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
            "thickness_calving_threshold": 50,
            "fractures": "false",
        },
    },
}

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
    "-d",
    "--distribution",
    dest="distribution",
    choices=dists.keys(),
    help="""Choose set.""",
    default="calving",
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
distribution_name = options.distribution

print(f"\nDrawing {n_samples} samples from distribution set {distribution_name}")
distributions = dists[distribution_name]["uq"]

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
    elif key == "frontal_melt_file":
        dist_sample[:, i] = [
            f"jib_ocean_forcing_{int(id)}_1980_2020.nc"
            for id in distributions[key].ppf(unif_sample[:, i])
        ]
    else:
        dist_sample[:, i] = distributions[key].ppf(unif_sample[:, i])


# Convert to Pandas dataframe, append column headers, output as csv
df = pd.DataFrame(dist_sample, columns=distributions.keys())
df.to_csv(outfile, index=True, index_label="id")

print("\nAdding default values\n")
for key, val in dists[distribution_name]["default_values"].items():
    if key not in df.columns:
        df[key] = val
        print(f"{key}: {val}")

df.to_csv(f"ensemble_{outfile}", index=True, index_label="id")
