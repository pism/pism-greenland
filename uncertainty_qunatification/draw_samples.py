#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import truncnorm, gamma, uniform, randint

from SALib.sample import saltelli
from SALib.analyze import sobol

parser = ArgumentParser()
parser.description = "Draw samples using the Saltelli methods"
parser.add_argument(
    "-s", "--n_samples", dest="n_samples", type=int, help="""number of samples to draw. default=10.""", default=10
)
parser.add_argument("OUTFILE", nargs=1, help="Ouput file (CSV)", default="velocity_calibration_samples.csv")
options = parser.parse_args()
n_samples = options.n_samples
outfile = options.OUTFILE[-1]


distributions = {
    "GCM": randint(0, 4),
    "FICE": truncnorm(-4 / 4.0, 4.0 / 4, loc=8, scale=4),
    "FSNOW": truncnorm(-4.1 / 3, 4.1 / 3, loc=4.1, scale=1.5),
    "PRS": uniform(loc=5, scale=2),
    "RFR": truncnorm(-0.4 / 0.3, 0.4 / 0.3, loc=0.5, scale=0.2),
    "OCM": randint(-1, 2),
    "OCS": randint(-1, 2),
    "TCT": randint(-1, 2),
    "VCM": truncnorm(-0.35 / 0.2, 0.35 / 0.2, loc=1, scale=0.2),
    "PPQ": truncnorm(-0.35 / 0.2, 0.35 / 0.2, loc=0.6, scale=0.2),
    "SIAE": gamma(1.5, scale=0.8, loc=1),
}

distributions = {
    "vcm": uniform(loc=0.5, scale=0.5),
    "fracture_softening": uniform(loc=0.01, scale=0.99),
    "fracture_rate": uniform(loc=0.1, scale=0.8),
    "fracture_threshold": uniform(loc=40e3, scale=110e3),
    "fracture_healing_rate": uniform(loc=0.0, scale=2.0),
    "fracture_healing_threshold": uniform(loc=1e-11, scale=9.9e-10),
}


# Generate the Sobol sequence samples with uniform distributions

problem = {
    "num_vars": len(distributions.keys()),
    "names": distributions.keys(),
    "bounds": [[0, 1]] * len(distributions.keys()),
}

# Generate samples
unif_sample = saltelli.sample(problem, n_samples, calc_second_order=False)


# To hold the transformed variables
dist_sample = np.zeros_like(unif_sample)


# For each variable, transform with the inverse of the CDF (inv(CDF)=ppf)
for i, key in enumerate(distributions.keys()):
    dist_sample[:, i] = distributions[key].ppf(unif_sample[:, i])


# Convert to Pandas dataframe, append column headers, output as csv
df = pd.DataFrame(dist_sample, columns=distributions.keys())

df["climate"] = "given"
df["hydrology"] = "routing"
df["frontal_melt"] = "discharge_routing"
df["climate_file"] = "DMI-HIRHAM5_GL2_ERAI_1980_2016_EPSG3413_4500M_DM.nc"
df["runoff_file"] = "DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_4500M_DM.nc"
df["frontal_melt_file"] = "jib_ocean_forcing_ctrl_1980_2020.nc"
df["calving_rate_scaling_file"] = "seasonal_calving.nc"
df["thickness_calving_threshold"] = 250
df["gamma_T"] = 1.25e-4
df["salinity"] = ""
df["pseudo_plastic_q"] = 0.6
df["sia_e"] = 1.25
df["fractures"] = "true"


df.to_csv(outfile, index=True, index_label="id")
