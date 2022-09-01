#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import uniform, randint, truncnorm, gamma

from SALib.sample import saltelli
from pyDOE import lhs

short2long = {
    "SIAE": "sia_e",
    "SSAN": "ssa_n",
    "PPQ": "pseudo_plastic_q",
    "TEFO": "till_effective_fraction_overburden",
    "PHIMIN": "phi_min",
    "PHIMAX": "phi_max",
    "ZMIN": "z_min",
    "ZMAX": "z_max",
}

pr2num = {
    0: "low-5",
    1: "low-mean",
    2: "low-95",
    3: "main-5",
    4: "main-mean",
    5: "main-95",
    6: "high-5",
    7: "high-mean",
    8: "high-95",
}

tas2num = {
    0: "main-5",
    1: "main-mean",
    2: "main-95",
    3: "S1-5",
    4: "S1-mean",
    5: "S1-95",
    6: "S2-5",
    7: "S2-mean",
    8: "S2-95",
    9: "S3-5",
    10: "S3-mean",
    11: "S3-95",
    12: "S4-5",
    13: "S4-mean",
    14: "S4-95",
}


dists = {
    "calving": {
        "uq": {
            "vcm": uniform(loc=250000, scale=250000),
        },
        "default_values": {
            "phi_min": 5,
            "phi_max": 40,
            "z_min": -700,
            "z_max": 700,
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.25,
            "till_effective_fraction_overburden": 0.02,
            "f_ice": 8,
            "f_snow": 3,
            "pr_paleo_file": "pr_Badgeley_etal_2020_id_main-mean.nc",
            "tas_paleo_file": "tas_Badgeley_etal_2020_id_main-mean.nc",
            "thickness_calving_threshold": 50,
            "eigen_calving_K": 1e19,
        },
    },
    "climate": {
        "uq": {
            "pr_paleo_file": randint(0, len(pr2num)),
            "tas_paleo_file": randint(0, len(tas2num)),
        },
        "default_values": {
            "phi_min": 5,
            "phi_max": 40,
            "z_min": -700,
            "z_max": 700,
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.25,
            "till_effective_fraction_overburden": 0.02,
            "f_ice": 8,
            "f_snow": 3,
            "thickness_calving_threshold": 50,
            "eigen_calving_K": 1e19,
            "vcm": 400000,
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
    "-m",
    "--method",
    dest="method",
    type=str,
    choices=["lhs", "saltelli"],
    help="""number of samples to draw. default=saltelli.""",
    default="saltelli",
)
parser.add_argument(
    "--posterior_file",
    help="Posterior predictive parameter file",
    default=None,
)
parser.add_argument(
    "OUTFILE",
    nargs=1,
    help="Ouput file (CSV)",
    default="velocity_calibration_samples.csv",
)
options = parser.parse_args()
n_draw_samples = options.n_samples
calc_second_order = options.calc_second_order
method = options.method
outfile = options.OUTFILE[-1]
distribution_name = options.distribution
posterior_file = options.posterior_file

print(f"\nDrawing {n_draw_samples} samples from distribution set {distribution_name}")
distributions = dists[distribution_name]["uq"]

problem = {
    "num_vars": len(distributions.keys()),
    "names": distributions.keys(),
    "bounds": [[0, 1]] * len(distributions.keys()),
}

keys_prior = list(distributions.keys())

# Generate uniform samples (i.e. one unit hypercube)
if method == "saltelli":
    unif_sample = saltelli.sample(
        problem, n_draw_samples, calc_second_order=calc_second_order
    )
elif method == "lhs":
    unif_sample = lhs(len(keys_prior), n_draw_samples)
else:
    print(f"Method {method} not available")

n_samples = unif_sample.shape[0]
# To hold the transformed variables
dist_sample = np.zeros_like(unif_sample, dtype=object)

# For each variable, transform with the inverse of the CDF (inv(CDF)=ppf)
for i, key in enumerate(keys_prior):
    if key == "pr_paleo_file":
        [
            print(f"pr_Badgeley_etal_2020_id_{pr2num[id]}.nc")
            for id in distributions[key].ppf(unif_sample[:, i])
        ]
        dist_sample[:, i] = [
            f"pr_Badgeley_etal_2020_id_{pr2num[id]}.nc"
            for id in distributions[key].ppf(unif_sample[:, i])
        ]
    elif key == "tas_paleo_file":
        dist_sample[:, i] = [
            f"tas_Badgeley_etal_2020_id_{tas2num[id]}.nc"
            for id in distributions[key].ppf(unif_sample[:, i])
        ]
    else:
        dist_sample[:, i] = distributions[key].ppf(unif_sample[:, i])

if posterior_file:
    X_posterior = pd.read_csv(posterior_file).drop(
        columns=["Unnamed: 0", "Model"], errors="ignore"
    )
    keys_mc = list(X_posterior.keys())
    keys = set(keys_prior + keys_mc)
    if len(keys_prior) + len(keys_mc) != len(keys):
        import sys

        print("Duplicate keys, exciting.")
        sys.exit()
    keys = keys_prior + keys_mc
    mc_indices = np.random.choice(range(X_posterior.shape[0]), n_samples)
    X_sample = X_posterior.to_numpy()[mc_indices, :]

    dist_sample = np.hstack((dist_sample, X_sample))

else:
    keys = keys_prior


# Convert to Pandas dataframe, append column headers, output as csv
df = pd.DataFrame(dist_sample, columns=keys)
df.to_csv(outfile, index=True, index_label="id")

print("\nAdding default values\n")
for key, val in dists[distribution_name]["default_values"].items():
    if key not in df.columns:
        df[key] = val
        print(f"{key}: {val}")

df.to_csv(f"ensemble_{outfile}", index=True, index_label="id")
