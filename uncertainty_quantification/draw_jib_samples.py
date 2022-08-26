#!/usr/bin/env python


from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scipy.stats.distributions import uniform, randint, truncnorm, gamma

from SALib.sample import saltelli


dists = {
    "fractures": {
        "uq": {
            "fracture_gamma": uniform(loc=0, scale=1),
            "fracture_gamma_h": uniform(loc=0, scale=1),
            "fracture_initiation_threshold": uniform(loc=40e3, scale=110e3),
            "healing_threshold": uniform(loc=1e-11, scale=9.9e-10),
            "fracture_softening": uniform(loc=0.65, scale=0.35),
            "vcm": uniform(loc=0.4, scale=0.4),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "frontal_melt_file": "jib_ocean_forcing_id_fjord_ctrl_1980_2020.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
            "gamma_T": 1.25e-4,
            "thickness_calving_threshold": 50,
            "fractures": "true",
            "parameter_a": 3e-4,
            "power_alpha": 0.39,
        },
    },
    "fractures_all": {
        "uq": {
            "fracture_gamma": uniform(loc=0, scale=1),
            "fracture_gamma_h": uniform(loc=0, scale=1),
            "fracture_initiation_threshold": uniform(loc=40e3, scale=110e3),
            "healing_threshold": uniform(loc=1e-11, scale=9.9e-10),
            "fracture_softening": uniform(loc=0.60, scale=0.40),
            "thickness_calving_threshold": uniform(loc=100, scale=400),
            "vcm": uniform(loc=0.4, scale=0.4),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "frontal_melt_file": "jib_ocean_forcing_id_fjord_ctrl_1980_2020.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
            "gamma_T": 1.25e-4,
            "fractures": "true",
            "parameter_a": 3e-4,
            "power_alpha": 0.39,
        },
    },
    "fractures_steady": {
        "uq": {
            "fracture_gamma": uniform(loc=0, scale=1),
            "fracture_gamma_h": uniform(loc=0, scale=1),
            "fracture_initiation_threshold": uniform(loc=40e3, scale=110e3),
            "healing_threshold": uniform(loc=1e-11, scale=9.9e-10),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "frontal_melt_file": "jib_ocean_forcing_ctrl_1980_2020.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
            "gamma_T": 1.25e-4,
            "thickness_calving_threshold": 300,
            "fractures": "true",
            "fracture_gamma": 0.4697265625,
            "fracture_gamma_h": 0.0,
            "fracture_softening": 0.8466796875,
            "fracture_initiation_threshold": 127548.828125,
            "healing_threshold": 4.3249023437500005e-10,
            "vcm": 0.5,
            "parameter_a": 3e-4,
            "power_alpha": 0.39,
        },
    },
    "fractures_melt": {
        "uq": {
            "fracture_softening": uniform(loc=0.90, scale=0.10),
            "gamma_T": uniform(loc=1.20e-4, scale=0.4e-4),
            "vcm": uniform(loc=0.45, scale=0.2),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "frontal_melt_file": "jib_ocean_forcing_id_fjord_ctrl_1980_2020.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
            "gamma_T": 1.25e-4,
            "thickness_calving_threshold": 300,
            "fractures": "true",
            "fracture_gamma": 0.4697265625,
            "fracture_gamma_h": 0.0,
            "fracture_softening": 0.8466796875,
            "fracture_initiation_threshold": 127548.828125,
            "healing_threshold": 4.3249023437500005e-10,
            "vcm": 0.5,
            "parameter_a": 3e-4,
            "power_alpha": 0.39,
        },
    },
    "calving": {
        "uq": {
            "vcm": uniform(loc=0.5, scale=0.2),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "frontal_melt_file": "jib_ocean_forcing_id_fjord_ctrl_1980_2020.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
            "gamma_T": 0.0001338671875,
            "thickness_calving_threshold": "thickness_calving_threshold_300_500_1998-01-01_2000-01-01.nc",
            "fractures": "true",
            "fracture_gamma": 0.4697265625,
            "fracture_gamma_h": 0.0,
            "fracture_softening": 0.9969726562500001,
            "fracture_initiation_threshold": 127548.828125,
            "healing_threshold": 4.3249023437500005e-10,
            "parameter_a": 3e-4,
            "power_alpha": 0.39,
        },
    },
    "ocean": {
        "uq": {
            "frontal_melt_file": randint(0, 10),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
            "gamma_T": 1.5e-4,
            "thickness_calving_threshold": 300,
            "fractures": "true",
            "fracture_gamma": 0.4697265625,
            "fracture_gamma_h": 0.0,
            "fracture_softening": 0.96767578125,
            "fracture_initiation_threshold": 127548.828125,
            "healing_threshold": 4.3249023437500005e-10,
            "vcm": 0.5,
            "parameter_a": 3e-4,
            "power_alpha": 0.39,
        },
    },
    "full": {
        "uq": {
            "frontal_melt_file": randint(0, 10),
            "vcm": uniform(loc=0.4, scale=0.6),
            "thickness_calving_threshold": uniform(loc=100, scale=400),
            "gamma_T": uniform(loc=1.00e-4, scale=0.6e-4),
            "power_alpha": uniform(loc=0.5, scale=0.5),
            "parameter_a": uniform(loc=1e-4, scale=5e-4),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "fractures": "false",
            "fracture_gamma": 0.4697265625,
            "fracture_gamma_h": 0.0,
            "fracture_softening": 0.96767578125,
            "fracture_initiation_threshold": 127548.828125,
            "healing_threshold": 4.3249023437500005e-10,
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "frontal_melt_file": "jib_ocean_forcing_id_fjord_ctrl_1980_2020.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
        },
    },
    "frontal_ablation": {
        "uq": {
            "frontal_melt_file": randint(0, 10),
            "vcm": uniform(loc=0.4, scale=0.6),
            "thickness_calving_threshold": uniform(loc=100, scale=400),
            "gamma_T": uniform(loc=1.00e-4, scale=0.6e-4),
            "power_alpha": uniform(loc=0.5, scale=0.5),
            "parameter_a": uniform(loc=1e-4, scale=5e-4),
            "calving_rate_scaling_file": randint(0, 2),
        },
        "default_values": {
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "fractures": "false",
            "fracture_gamma": 0.4697265625,
            "fracture_gamma_h": 0.0,
            "fracture_softening": 0.96767578125,
            "fracture_initiation_threshold": 127548.828125,
            "healing_threshold": 4.3249023437500005e-10,
            "climate": "given",
            "hydrology": "routing",
            "frontal_melt": "discharge_routing",
            "frontal_melt_file": "jib_ocean_forcing_id_fjord_ctrl_1980_2020.nc",
            "climate_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "runoff_file": "DMI-HIRHAM5_ERA_1980_2020_EPSG3413_4500M_DM.nc",
            "salinity": "",
            "pseudo_plastic_q": 0.6,
            "sia_e": 1.25,
            "ssa_n": 3.0,
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
    default="all",
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
