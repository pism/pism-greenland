#!/usr/bin/env python
# Copyright (C) 2016-21 Andy Aschwanden

import itertools
from collections import OrderedDict
import csv
import numpy as np
import os
import sys
import shlex
from os.path import join, abspath, realpath, dirname
import pandas as pd

try:
    import subprocess32 as sub
except:
    import subprocess as sub

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys


def current_script_directory():
    import inspect

    filename = inspect.stack(0)[0][1]
    return realpath(dirname(filename))


script_directory = current_script_directory()

sys.path.append(join(script_directory, "../resources"))
from resources import *


def map_dict(val, mdict):
    try:
        return mdict[val]
    except:
        return val


grid_choices = [
    18000,
    9000,
    6000,
    4500,
    3600,
    3000,
    2400,
    1800,
    1500,
    1200,
    900,
    600,
    450,
    300,
    150,
]

# set up the option parser
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Generating scripts for warming experiments."
parser.add_argument("FILE", nargs=1, help="Input file to restart from", default=None)
parser.add_argument(
    "-n",
    "--n_procs",
    dest="n",
    type=int,
    help="""number of cores/processors. default=140.""",
    default=140,
)
parser.add_argument(
    "-w",
    "--wall_time",
    dest="walltime",
    help="""walltime. default: 100:00:00.""",
    default="100:00:00",
)
parser.add_argument(
    "-q",
    "--queue",
    dest="queue",
    choices=list_queues(),
    help="""queue. default=long.""",
    default="long",
)
parser.add_argument(
    "-d",
    "--domain",
    dest="domain",
    choices=["gris", "gris_ext", "jib", "qaanaaq"],
    help="sets the modeling domain",
    default="gris",
)
parser.add_argument(
    "--exstep",
    dest="exstep",
    help="Writing interval for spatial time series",
    default=1,
)
parser.add_argument(
    "-f",
    "--o_format",
    dest="oformat",
    choices=["netcdf3", "netcdf4_parallel", "netcdf4_serial", "pnetcdf"],
    help="output format",
    default="netcdf4_parallel",
)
parser.add_argument(
    "-g",
    "--grid",
    dest="grid",
    type=int,
    choices=grid_choices,
    help="horizontal grid resolution",
    default=1500,
)
parser.add_argument(
    "--i_dir",
    dest="input_dir",
    help="input directory",
    default=abspath(join(script_directory, "..")),
)
parser.add_argument(
    "--o_dir", dest="output_dir", help="output directory", default="test_dir"
)
parser.add_argument(
    "--o_size",
    dest="osize",
    choices=["small", "medium", "big", "big_2d", "custom"],
    help="output size type",
    default="custom",
)
parser.add_argument(
    "-s",
    "--system",
    dest="system",
    choices=list_systems(),
    help="computer system to use.",
    default="pleiades_broadwell",
)
parser.add_argument(
    "-b",
    "--bed_type",
    dest="bed_type",
    choices=list_bed_types(),
    help="output size type",
    default="wc",
)
parser.add_argument(
    "--spatial_ts",
    dest="spatial_ts",
    choices=["basic", "standard", "none", "hydro", "calib"],
    help="output size type",
    default="calib",
)
parser.add_argument(
    "--hydrology",
    dest="hydrology",
    choices=["null", "diffuse", "routing", "distributed"],
    help="Basal hydrology model.",
    default="diffuse",
)
parser.add_argument(
    "-p",
    "--params",
    dest="params_list",
    help="Comma-separated list with params for sensitivity",
    default=None,
)
parser.add_argument(
    "--stable_gl",
    dest="float_kill_calve_near_grounding_line",
    action="store_false",
    help="Stable grounding line",
    default=True,
)
parser.add_argument(
    "--hot_spot",
    dest="hot_spot",
    action="store_true",
    help="Use hotpsot",
    default=False,
)
parser.add_argument(
    "--stress_balance",
    dest="stress_balance",
    choices=["sia", "ssa+sia", "ssa", "blatter"],
    help="stress balance solver",
    default="ssa+sia",
)
parser.add_argument(
    "--dataset_version",
    dest="version",
    choices=["2", "3", "3a", "4", "1980", "1980v3", "5_RAGIS", "2022"],
    help="Input data set version",
    default="2022",
)
parser.add_argument(
    "--vertical_velocity_approximation",
    dest="vertical_velocity_approximation",
    choices=["centered", "upstream"],
    help="How to approximate vertical velocities",
    default="upstream",
)
parser.add_argument(
    "-e",
    "--ensemble_file",
    dest="ensemble_file",
    help="File that has all combinations for ensemble study",
    default="../uncertainty_quantification/initialization.csv",
)
parser.add_argument(
    "-L",
    "--comp_level",
    dest="compression_level",
    help="Compression level for output file. Only works with netcdf4_serial.",
    default=2,
)

parser.add_argument(
    "--start_year", dest="start_year", type=int, help="Simulation start year", default=0
)
parser.add_argument(
    "--duration", dest="duration", type=int, help="Years to simulate", default=50
)
parser.add_argument(
    "--step", dest="step", type=int, help="Step in years for restarting", default=50
)

options = parser.parse_args()

nn = options.n
input_dir = abspath(options.input_dir)
output_dir = abspath(options.output_dir)
spatial_tmp_dir = abspath(options.output_dir + "_tmp")

compression_level = options.compression_level
oformat = options.oformat
osize = options.osize
queue = options.queue
walltime = options.walltime
system = options.system

ensemble_file = options.ensemble_file

spatial_ts = options.spatial_ts

bed_type = options.bed_type
climate = "flux"
exstep = options.exstep
float_kill_calve_near_grounding_line = options.float_kill_calve_near_grounding_line
grid = options.grid
hydrology = options.hydrology
stress_balance = options.stress_balance
vertical_velocity_approximation = options.vertical_velocity_approximation
version = options.version
ocean = "const"
hot_spot = options.hot_spot

domain = options.domain
pism_exec = generate_domain(domain)

if options.FILE is None:
    print("Missing input file")
    import sys

    sys.exit()
else:
    input_file = options.FILE[0]

if domain.lower() in ("greenland_ext", "gris_ext"):
    pism_dataname = (
        "$input_dir/data_sets/bed_dem/pism_Greenland_ext_{}m_mcb_jpl_v{}_{}.nc".format(
            grid, version, bed_type
        )
    )
else:
    pism_dataname = (
        "$input_dir/data_sets/bed_dem/pism_Greenland_{}m_mcb_jpl_v{}_{}.nc".format(
            grid, version, bed_type
        )
    )

climate_file = "$input_dir/data_sets/climate_forcing/DMI-HIRHAM5_GL2_ERAI_2001_2014_YDM_BIL_EPSG3413_{}m.nc".format(
    grid
)

# regridvars = "litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk"
regridvars = "litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume"

dirs = {"output": "$output_dir", "spatial_tmp": "$spatial_tmp_dir"}
for d in ["performance", "state", "scalar", "spatial", "jobs", "basins"]:
    dirs[d] = "$output_dir/{dir}".format(dir=d)

if spatial_ts == "none":
    del dirs["spatial"]

# use the actual path of the run scripts directory (we need it now and
# not during the simulation)
scripts_dir = join(output_dir, "run_scripts")
if not os.path.isdir(scripts_dir):
    os.makedirs(scripts_dir)

# use the actual path of the experiment table directory (we need it now and
# not during the simulation)
exp_dir = join(output_dir, "experiments")
if not os.path.isdir(exp_dir):
    os.makedirs(exp_dir)

# generate the config file *after* creating the output directory
pism_config = "calibrate_config"
pism_config_nc = join(output_dir, pism_config + ".nc")

cmd = "ncgen -o {output} {input_dir}/config/{config}.cdl".format(
    output=pism_config_nc, input_dir=input_dir, config=pism_config
)
sub.call(shlex.split(cmd))

# these Bash commands are added to the beginning of the run scrips
run_header = """# stop if a variable is not defined
set -u
# stop on errors
set -e

# path to the config file
config="{config}"
# path to the input directory (input data sets are contained in this directory)
input_dir="{input_dir}"
# output directory
output_dir="{output_dir}"
# temporary directory for spatial files
spatial_tmp_dir="{spatial_tmp_dir}"

# create required output directories
for each in {dirs};
do
  mkdir -p $each
done

""".format(
    input_dir=input_dir,
    output_dir=output_dir,
    spatial_tmp_dir=spatial_tmp_dir,
    config=pism_config_nc,
    dirs=" ".join(list(dirs.values())),
)
if system not in ("debug", "pleiades"):
    cmd = "lfs setstripe -c -1 {}".format(dirs["output"])
    sub.call(shlex.split(cmd))
    cmd = "lfs setstripe -c -1 {}".format(dirs["spatial_tmp"])
    sub.call(shlex.split(cmd))

# ########################################################
# set up parameters
# ########################################################

ssa_e = 1.0
tlftw = 0.1

uq_df = pd.read_csv(ensemble_file)
uq_df.fillna(False, inplace=True)


tsstep = "yearly"

scripts = []
scripts_combinded = []

simulation_start_year = options.start_year
simulation_end_year = options.start_year + options.duration
restart_step = options.step

if restart_step > (simulation_end_year - simulation_start_year):
    print("Error:")
    print(
        (
            "restart_step > (simulation_end_year - simulation_start_year): {} > {}".format(
                restart_step, simulation_end_year - simulation_start_year
            )
        )
    )
    print("Try again")
    import sys

    sys.exit(0)

batch_header, batch_system = make_batch_header(system, nn, walltime, queue)

for n, row in enumerate(uq_df.iterrows()):
    combination = row[1]
    print(combination)
    phi_min = combination["phi_min"]
    phi_max = combination["phi_max"]
    z_min = combination["z_min"]
    z_max = combination["z_max"]
    ttphi = "{},{},{},{}".format(phi_min, phi_max, z_min, z_max)

    vversion = "v" + str(version)

    name_options = OrderedDict()
    m_id = combination["id"]
    try:
        m_id = int(m_id)
    except:
        pass
    name_options["id"] = m_id

    full_exp_name = "_".join(
        [
            vversion,
            "_".join(["_".join([k, str(v)]) for k, v in list(name_options.items())]),
        ]
    )
    full_outfile = "g{grid}m_{experiment}.nc".format(
        grid=grid, experiment=full_exp_name
    )

    # All runs in one script file for coarse grids that fit into max walltime
    script_combined = join(scripts_dir, "cc_g{}m_{}_j.sh".format(grid, full_exp_name))
    with open(script_combined, "w") as f_combined:

        outfiles = []
        job_no = 0
        for start in range(simulation_start_year, simulation_end_year, restart_step):
            job_no += 1

            end = start + restart_step

            experiment = "_".join(
                [
                    vversion,
                    "_".join(
                        ["_".join([k, str(v)]) for k, v in list(name_options.items())]
                    ),
                    "{}".format(start),
                    "{}".format(end),
                ]
            )

            script = join(scripts_dir, "cc_g{}m_{}.sh".format(grid, experiment))
            scripts.append(script)

            for filename in script:
                try:
                    os.remove(filename)
                except OSError:
                    pass

            if start == simulation_start_year:
                f_combined.write(batch_header)
                f_combined.write(run_header)

            with open(script, "w") as f:

                f.write(batch_header)
                f.write(run_header)

                outfile = "{domain}_g{grid}m_{experiment}.nc".format(
                    domain=domain.lower(), grid=grid, experiment=experiment
                )

                pism = generate_prefix_str(pism_exec)

                general_params_dict = {
                    "profile": join(
                        dirs["performance"],
                        "profile_${job_id}.py".format(**batch_system),
                    ),
                    "ys": start,
                    "ye": end,
                    "calendar": "365_day",
                    "o": join(dirs["state"], outfile),
                    "o_format": oformat,
                    "output.compression_level": compression_level,
                    "config_override": "$config",
                    "energy.ch_warming.enabled": "false",
                }

                if start == simulation_start_year:
                    general_params_dict["bootstrap"] = ""
                    general_params_dict["i"] = pism_dataname
                    if hot_spot:
                        general_params_dict[
                            "energy.bedrock_thermal.file"
                        ] = f"$input_dir/data_sets/bheatflux/Geothermal_Heat_Flux_Greenland_corrected_g{grid}m.nc"
                    general_params_dict["regrid_file"] = input_file
                    general_params_dict["regrid_vars"] = regridvars
                else:
                    general_params_dict["i"] = regridfile

                if osize != "custom":
                    general_params_dict["o_size"] = osize
                else:
                    general_params_dict[
                        "output.sizes.medium"
                    ] = "sftgif,velsurf_mag,tempicethk_basal,velsurf,velbase_mag"

                if start == simulation_start_year:
                    grid_params_dict = generate_grid_description(grid, domain)
                else:
                    grid_params_dict = generate_grid_description(
                        grid, domain, restart=True
                    )

                sliding_law = "pseudo_plastic"
                if "sliding_law" in combination:
                    sliding_law = combination["sliding_law"]

                u_threshold = 100.0
                if "u_threshold" in combination:
                    u_threshold = combination["u_threshold"]
                sb_params_dict = {
                    sliding_law: "",
                    "sia_e": combination["sia_e"],
                    "ssa_n": combination["ssa_n"],
                    "ssa_e": 1.0,
                    "pseudo_plastic_q": combination["pseudo_plastic_q"],
                    "till_effective_fraction_overburden": combination[
                        "till_effective_fraction_overburden"
                    ],
                    "vertical_velocity_approximation": vertical_velocity_approximation,
                    "stress_balance.blatter.enhancement_factor": combination["sia_e"],
                    "basal_yield_stress.mohr_coulomb.till_log_factor_transportable_water": tlftw,
                    "basal_resistance.pseudo_plastic.u_threshold": u_threshold,
                }

                if start == simulation_start_year:
                    sb_params_dict["topg_to_phi"] = ttphi

                # If stress balance choice is made in file, overwrite command line option
                stress_balance_params_dict = generate_stress_balance(
                    stress_balance, sb_params_dict
                )

                climate_params_dict = generate_climate(
                    climate, force_to_thickness_file=pism_dataname
                )

                hydro_params_dict = generate_hydrology(hydrology)

                calving_params_dict = generate_calving(
                    "vonmises_calving", front_retreat_file=pism_dataname
                )
                scalar_ts_dict = generate_scalar_ts(
                    outfile,
                    tsstep,
                    start=simulation_start_year,
                    end=simulation_end_year,
                    odir=dirs["scalar"],
                )

                all_params_dict = merge_dicts(
                    general_params_dict,
                    grid_params_dict,
                    stress_balance_params_dict,
                    climate_params_dict,
                    hydro_params_dict,
                    calving_params_dict,
                    scalar_ts_dict,
                )

                if not spatial_ts == "none":
                    exvars = spatial_ts_vars[spatial_ts]
                    spatial_ts_dict = generate_spatial_ts(
                        outfile, exvars, exstep, odir=dirs["spatial_tmp"], split=False
                    )

                    all_params_dict = merge_dicts(all_params_dict, spatial_ts_dict)
                if stress_balance == "blatter":
                    del all_params_dict["skip"]
                    all_params_dict["time_stepping.adaptive_ratio"] = 25

                all_params = " \\\n  ".join(
                    ["-{} {}".format(k, v) for k, v in list(all_params_dict.items())]
                )

                if system == "debug":
                    redirect = " 2>&1 | tee {jobs}/job_{job_no}.${job_id}"
                else:
                    redirect = " > {jobs}/job_{job_no}.${job_id} 2>&1"

                template = "{mpido} {pism} {params}" + redirect

                context = merge_dicts(
                    batch_system,
                    dirs,
                    {"job_no": job_no, "pism": pism, "params": all_params},
                )
                cmd = template.format(**context)

                f.write(cmd)
                f.write("\n")
                f.write(batch_system.get("footer", ""))

                f_combined.write(cmd)
                f_combined.write("\n\n")

                f_combined.write("\n")
                f_combined.write("\n")
                if not spatial_ts == "none":
                    f_combined.write(
                        "mv {tmpfile} {ofile}\n".format(
                            tmpfile=spatial_ts_dict["extra_file"],
                            ofile=join(dirs["spatial"], "ex_" + outfile),
                        )
                    )
                    f_combined.write("\n")

                regridfile = join(dirs["state"], outfile)
                outfiles.append(outfile)

            f_combined.write(batch_system.get("footer", ""))

        scripts_combinded.append(script_combined)

scripts = uniquify_list(scripts)
scripts_combinded = uniquify_list(scripts_combinded)
print("\n".join([script for script in scripts]))
print("\nwritten\n")
print("\n".join([script for script in scripts_combinded]))
print("\nwritten\n")
