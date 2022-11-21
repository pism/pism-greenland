#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import uniform, randint, truncnorm, gamma

from SALib.sample import saltelli
from pyDOE2 import lhs

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

dists = {
    "gris_pseudo_plastic": {
        "uq": {
            "sia_e": uniform(
                loc=1.0, scale=3.0
            ),  # uniform between 1 and 4    AS16 best value: 1.25
            "sia_n": uniform(
                loc=3.0, scale=0.5
            ),  # uniform between 3 and 3.5  AS16 best value: 3.25
            "ssa_n": uniform(
                loc=3.0, scale=0.5
            ),  # uniform between 3 and 3.5  AS16 best value: 3.25
            "pseudo_plastic_q": uniform(
                loc=0.25, scale=0.7
            ),  # uniform between 0.25 and 0.95
            "pseudo_plastic_uthreshold": uniform(loc=50.0, scale=150.0),
            "till_effective_fraction_overburden": uniform(
                loc=0.01, scale=0.02
            ),  # uniform between 0.015 and 0.030
            "phi_min": uniform(loc=5.0, scale=15.0),  # uniform between  5 and 20
            "z_min": uniform(loc=-1000, scale=1000),  # uniform between -1000 and 0
            "z_max": uniform(loc=0, scale=1000),  # uniform between 0 and 1000
        },
        "default_values": {
            "phi_max": 40,
            "sliding_law": "pseudo_plastic",
        },
    },
    "gris_regularized_coulomb": {
        "uq": {
            "sia_e": uniform(
                loc=1.0, scale=3.0
            ),  # uniform between 1 and 4    AS16 best value: 1.25
            "sia_n": uniform(
                loc=3.0, scale=1.0
            ),  # uniform between 3 and 3.5  AS16 best value: 3.25
            "ssa_n": uniform(
                loc=3.0, scale=0.5
            ),  # uniform between 3 and 3.5  AS16 best value: 3.25
            "pseudo_plastic_q": uniform(
                loc=0.25, scale=0.7
            ),  # uniform between 0.25 and 0.95
            "pseudo_plastic_uthreshold": uniform(loc=50.0, scale=150.0),
            "till_effective_fraction_overburden": uniform(
                loc=0.01, scale=0.02
            ),  # uniform between 0.015 and 0.030
            "phi_min": uniform(loc=5.0, scale=15.0),  # uniform between  5 and 20
            "z_min": uniform(loc=-1000, scale=1000),  # uniform between -1000 and 0
            "z_max": uniform(loc=0, scale=1000),  # uniform between 0 and 1000
        },
        "default_values": {
            "phi_max": 40,
            "sliding_law": "regularized_coulomb",
        },
    },
}

parser = ArgumentParser()
parser.description = "Generate UQ using Latin Hypercube or Sobol Sequences."
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
    default="gris_pseudo_plastic",
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
    default="lhs",
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
    unif_sample = lhs(len(keys_prior), n_draw_samples, random_state=2)
else:
    print(f"Method {method} not available")

n_samples = unif_sample.shape[0]
# To hold the transformed variables
dist_sample = np.zeros_like(unif_sample, dtype="object")

sb_dict = {0: "ssa+sia", 1: "blatter"}
# For each variable, transform with the inverse of the CDF (inv(CDF)=ppf)
for i, key in enumerate(keys_prior):
    if key == "stress_balance":
        dist_sample[:, i] = [
            f"{sb_dict[int(id)]}" for id in distributions[key].ppf(unif_sample[:, i])
        ]
    else:
        dist_sample[:, i] = distributions[key].ppf(unif_sample[:, i])


# Convert to Pandas dataframe, append column headers, output as csv
df = pd.DataFrame(dist_sample, columns=keys_prior)
df.to_csv(outfile, index=True, index_label="id")

print("\nAdding default values\n")
for key, val in dists[distribution_name]["default_values"].items():
    if key not in df.columns:
        df[key] = val
        print(f"{key}: {val}")

df.to_csv(f"ensemble_{outfile}", index=True, index_label="id")
