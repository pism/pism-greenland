#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import uniform, randint, truncnorm, gamma

from SALib.sample import saltelli
from pyDOE2 import lhs


dists = {
    "ctrl": {
        "uq": {
        },
        "default_values": {
            "sia_e": 2.608046,
            "sia_n": 3.0,
            "ssa_n": 3.309718,
            "pseudo_plastic_q": 0.7508221,
            "pseudo_plastic_uthreshold": 100,
            "till_effective_fraction_overburden": 0.01845403,
            "phi_min": 7.193718,
            "z_min": -369.6359,
            "z_max": 243.8239,
            "phi_max": 42.79528,
            "sliding_law": "pseudo_plastic",
            "climate": "paleo_searise",
            "calving.eigen_calving.K": 2.5e18,
            "f_ice": 8,
            "f_snow": 3,
            "refreeze": 0.6,
            "std_dev": 5,
            "pr_paleo_file": "pism_dT.nc",
            "tas_paleo_file": "pism_dT.nc",
            "calving.thickness_calving.threshold": 100,
            "vcm": 0.5,
            "atmosphere.elevation_change.temperature_lapse_rate": 6
        },
    },
    "climate": {
        "uq": { "f_ice": uniform(6, 6),
                "f_snow": uniform(2, 3),
                "std_dev": uniform(3, 4),
                "refreeze": uniform(0.4, 0.4),
                "atmosphere.precip_exponential_factor_for_temperature": uniform(0.05, 0.03),
                "sia_n": uniform(1, 3),                
        },
        "default_values": {
            "sia_e": 2.608046,
            "ssa_n": 3.309718,
            "pseudo_plastic_q": 0.7508221,
            "pseudo_plastic_uthreshold": 100,
            "till_effective_fraction_overburden": 0.01845403,
            "phi_min": 7.193718,
            "z_min": -369.6359,
            "z_max": 243.8239,
            "phi_max": 42.79528,
            "sliding_law": "pseudo_plastic",
            "calving.eigen_calving.K": 2.5e18,
            "f_ice": 8,
            "f_snow": 3,
            "pr_paleo_file": "pism_dT.nc",
            "tas_paleo_file": "pism_dT.nc",
            "calving.thickness_calving.threshold": 100,
            "vcm": 0.5,
        },
    },
    "climate-calving": {
        "uq": { "f_ice": uniform(6, 6),
                "f_snow": uniform(2, 3),
                "std_dev": uniform(3, 4),
                "refreeze": uniform(0.4, 0.4),
                "vcm": uniform(0.4, 0.4),
                "calving.eigen_calving.K": uniform(1e16, 9.99e18),
                "sia_n": uniform(1, 3),
        },
        "default_values": {
            "sia_e": 2.608046,
            "ssa_n": 3.309718,
            "pseudo_plastic_q": 0.7508221,
            "pseudo_plastic_uthreshold": 100,
            "till_effective_fraction_overburden": 0.01845403,
            "phi_min": 7.193718,
            "z_min": -369.6359,
            "z_max": 243.8239,
            "phi_max": 42.79528,
            "sliding_law": "pseudo_plastic",
            "climate": "paleo_searise",
            "calving.eigen_calving.K": 2.5e18,
            "f_ice": 8,
            "f_snow": 3,
            "pr_paleo_file": "pism_dT.nc",
            "tas_paleo_file": "pism_dT.nc",
            "calving.thickness_calving.threshold": 100,
            "vcm": 0.5,
            "atmosphere.elevation_change.temperature_lapse_rate": 6
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
    default="climate",
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
    help="""number of samples to draw. default=lhs.""",
    default="lhs",
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
