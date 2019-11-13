#!/usr/bin/env python

from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import uniform, randint
from collections import OrderedDict

from SALib.sample import saltelli
from SALib.analyze import sobol

parser = ArgumentParser()
parser.description = "Draw samples using the Saltelli methods"
parser.add_argument(
    "-s", "--n_samples", dest="n_samples", type=int, help="""number of samples to draw. default=10.""", default=10
)
parser.add_argument("OUTFILE", nargs=1, help="Ouput file (CSV)", default="saltelli_samples.csv")
options = parser.parse_args()
n_samples = options.n_samples
outfile = options.OUTFILE[-1]

distributions = OrderedDict()
distributions["m_min"] = uniform(loc=-1.5, scale=0.5)
distributions["m_max"] = uniform(loc=4.0, scale=1.0)
distributions["h_min"] = randint(50, 150)
distributions["h_ela"] = randint(1500, 1800)
distributions["h_max"] = randint(2500, 3000)

# Names of all the variables
keys = distributions.keys()
print(distributions)
# Generate the Sobol sequence samples with uniform distributions

problem = {"num_vars": len(keys), "names": keys, "bounds": [[0, 1]] * len(keys)}

# Generate samples
unif_sample = saltelli.sample(problem, n_samples, calc_second_order=False)


# To hold the transformed variables
dist_sample = np.zeros_like(unif_sample)


# For each variable, transform with the inverse of the CDF (inv(CDF)=ppf)
for i, key in enumerate(keys):
    dist_sample[:, i] = distributions[key].ppf(unif_sample[:, i])


# Save to CSV file using Pandas DataFrame and to_csv method
header = keys
# Convert to Pandas dataframe, append column headers, output as csv
df = pd.DataFrame(data=dist_sample, columns=header)
df.to_csv(outfile, index=True, index_label="id")
