#!/usr/bin/env python
# Copyright (C) 2019-24 Andy Aschwanden

# Initialize Greenland LGM

import itertools
from collections import OrderedDict
import numpy as np
import os
import sys
import shlex
from os.path import join, abspath, realpath, dirname
import pandas as pd
import xarray as xr
from typing import Any, Dict, List, Union

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
    15000,
    12000,
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
    1000,
]

# set up the option parser
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Generating scripts for warming experiments."
parser.add_argument(
    "-n",
    "--n_procs",
    dest="n",
    type=int,
    help="""number of cores/processors. default=48.""",
    default=48,
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
    default="t2small",
)
parser.add_argument(
    "--options",
    dest="commandline_options",
    help="""Here you can add command-line options""",
)
parser.add_argument(
    "-d",
    "--domain",
    dest="domain",
    choices=["gris", "gris_ext"],
    help="sets the modeling domain",
    default="gris_ext",
)
parser.add_argument(
    "--exstep",
    dest="exstep",
    help="Writing interval for spatial time series",
    default=100,
)
parser.add_argument(
    "--tsstep",
    dest="tsstep",
    help="Writing interval for scalar time series",
    default="yearly",
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
    "-L",
    "--compression_level",
    dest="compression_level",
    type=int,
    help="Compression level for output file.",
    default=2,
)
parser.add_argument(
    "-g",
    "--grid",
    dest="grid",
    type=int,
    choices=grid_choices,
    help="horizontal grid resolution",
    default=4500,
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
    "-i",
    "--initial_state_file",
    dest="initialstatefile",
    help="Input file to restart from",
    default=None,
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
    choices=["basic", "paleo", "paleo_tracer"],
    help="output size type",
    default="paleo",
)
parser.add_argument(
    "--hydrology",
    dest="hydrology",
    choices=["routing", "routing_steady", "diffuse"],
    help="Basal hydrology model.",
    default="diffuse",
)
parser.add_argument(
    "--calving",
    dest="calving",
    choices=[
        "float_kill",
        "vonmises_calving",
        "eigen_calving",
        "hayhurst_calving",
        "hybrid_calving",
        "thickness_calving",
    ],
    help="Choose calving law",
    default="hybrid_calving",
)
parser.add_argument(
    "--stable_gl",
    dest="float_kill_calve_near_grounding_line",
    action="store_false",
    help="Stable grounding line",
    default=True,
)
parser.add_argument(
    "--gid",
    dest="gid",
    choices=["s2457", "s2524"],
    help="Choose GID for Pleiades",
    default="s2524",
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
    choices=["2023_GRIMP"],
    help="input data set version",
    default="2023_GRIMP",
)
parser.add_argument(
    "--bed_def",
    dest="bed_def",
    choices=["off", "lc"],
    help="Bedrock deformation model. Default=lc",
    default="lc",
)
parser.add_argument(
    "--age",
    action="store_true",
    help="Calculate age field. Default=False",
    default=False,
)
parser.add_argument("--start", help="Simulation start year", default=-125000)
parser.add_argument("--end", help="Simulation end year", default=0)
parser.add_argument(
    "-e",
    "--ensemble_file",
    dest="ensemble_file",
    help="File that has all combinations for ensemble study",
    default=None,
)

parser.add_argument(
    "--data_dir",
    dest="data_dir",
    help="data directory",
    default=abspath(join(script_directory, "../data_sets/")),
)

options = parser.parse_args()
commandline_options = options.commandline_options

start_date = options.start
end_date = options.end

nn = options.n
input_dir = abspath(options.input_dir)
data_dir = abspath(options.data_dir)
output_dir = abspath(options.output_dir)

compression_level = options.compression_level
oformat = options.oformat
osize = options.osize
queue = options.queue
walltime = options.walltime
system = options.system
gid = options.gid

age = options.age
initialstatefile = options.initialstatefile
spatial_ts = options.spatial_ts
bed_type = options.bed_type
exstep = options.exstep
tsstep = options.tsstep
float_kill_calve_near_grounding_line = options.float_kill_calve_near_grounding_line
grid = grid_resolution = options.grid
hydrology = options.hydrology
bed_def = options.bed_def
stress_balance = options.stress_balance
version = options.version

ensemble_file = options.ensemble_file

domain = options.domain
pism_exec = generate_domain(domain)

pism_dataname = False
if domain.lower() in ("greenland_ext", "gris_ext"):
    pism_dataname = (
        f"$data_dir/bed_dem/pism_Greenland_ext_{grid}m_v{version}_{bed_type}.nc"
    )
else:
    pism_dataname = (
        f"$data_dir/bed_dem/pism_Greenland_{grid}m_v{version}_{bed_type}.nc"
    )

#regridvars = "litho_temp,enthalpy,tillwat,bmelt,ice_area_specific_volume,thk,topg,isochrone_depth,isochronal_layer_thickness"
regridvars = "litho_temp,enthalpy,tillwat,bmelt,ice_area_specific_volume,thk,topg"

master_config_file = get_path_to_config()


dirs = {"output": "$output_dir"}
for d in ["performance", "state", "scalar", "snap", "spatial", "jobs", "basins"]:
    dirs[d] = f"$output_dir/{d}"

if spatial_ts == "none":
    del dirs["spatial"]

# use the actual path of the run scripts directory (we need it now and
# not during the simulation)
scripts_dir = join(output_dir, "run_scripts")
if not os.path.isdir(scripts_dir):
    os.makedirs(scripts_dir)

# use the actual path of the time file directory (we need it now and
# not during the simulation)
time_dir = join(output_dir, "time_forcing")
if not os.path.isdir(time_dir):
    os.makedirs(time_dir)

# use the actual path of the uq directory
uq_dir = join(output_dir, "uq")
if not os.path.isdir(uq_dir):
    os.makedirs(uq_dir)

# generate the config file *after* creating the output directory
pism_config = "paleo"
pism_config_nc = join(output_dir, pism_config + ".nc")

cmd = f"ncgen -o {pism_config_nc} {input_dir}/config/{pism_config}.cdl"
sub.call(shlex.split(cmd))

m_dirs = " ".join(list(dirs.values()))

# these Bash commands are added to the beginning of the run scrips
run_header = f"""# stop if a variable is not defined
set -u
# stop on errors
set -e

# path to the config file
config="{pism_config_nc}"
# path to the input directory (input data sets are contained in this directory)
input_dir="{input_dir}"
# path to data directory
data_dir="{data_dir}"
# output directory
output_dir="{output_dir}"

# create required output directories
for each in {m_dirs};
    do
      mkdir -p $each
done\n\n
"""


ensemble_infile = os.path.split(ensemble_file)[-1]
ensemble_outfile = join(uq_dir, ensemble_infile)

cmd = f"cp {ensemble_file} {ensemble_outfile}"
sub.call(shlex.split(cmd))

ensemble_infile_2 = ensemble_infile.split("ensemble_")[-1]
ensemble_outfile_2 = join(uq_dir, ensemble_infile_2)

cmd = f"cp {ensemble_file} {ensemble_outfile}"
sub.call(shlex.split(cmd))


# ########################################################
# set up model initialization
# ########################################################

ssa_n = 3.25
ssa_e = 1.0

uq_df = pd.read_csv(ensemble_file)
uq_df.fillna(False, inplace=True)

scripts = []
scripts_post = []

simulation_start_year = options.start
simulation_end_year = options.end

batch_header, batch_system = make_batch_header(system, nn, walltime, queue, gid=gid)
post_header = make_batch_post_header(system)

for n, row in enumerate(uq_df.iterrows()):
    combination = row[1]
    print(combination)

    name_options = {}
    try:
        name_options["id"] = int(combination["id"])
    except:
        name_options["id"] = combination["id"]

    vversion = "v" + str(version)
    full_exp_name = "_".join(
        [
            vversion,
            "_".join(["_".join([k, str(v)]) for k, v in list(name_options.items())]),
        ]
    )
    full_outfile = "g{grid}m_{experiment}.nc".format(
        grid=grid, experiment=full_exp_name
    )

    experiment = "_".join(
        [
            vversion,
            "_".join(["_".join([k, str(v)]) for k, v in list(name_options.items())]),
            "{}".format(start_date),
            "{}".format(end_date),
        ]
    )

    script = join(scripts_dir, f"{domain}_g{grid_resolution}m_{experiment}.sh")
    scripts.append(script)

    for filename in script:
        try:
            os.remove(filename)
        except OSError:
            pass

    with open(script, "w") as f:

        f.write(batch_header)
        f.write(run_header)

        pism = generate_prefix_str(pism_exec)

        general_params_dict = {
            "profile": join(
                dirs["performance"], "profile_${job_id}.py".format(**batch_system)
            ),
            "time.start": start_date,
            "time.end": end_date,
            "time.calendar": "365_day",
            "input.forcing.time_extrapolation": "true",
            "energy.bedrock_thermal.file": "$data_dir/bheatflux/Geothermal_heatflux_map_v2.1_g450m.nc",
            "output.format": oformat,
            "output.compression_level": compression_level,
            "config_override": "$config",
            "stress_balance.ice_free_thickness_standard": 5,
        }
        
        if bed_def != "off":
            general_params_dict["bed_deformation.model"] = bed_def

        if age:
            equally_spaced_layers = -np.arange(start_date + 5000, end_date, 5000)
            select_layers = -start_date - np.array([3, 8, 9, 11.7, 12.8, 14.7, 19, 29, 57]) * 1_000
            all_layers = " ".join([str(x) for x in np.sort(np.hstack([equally_spaced_layers, select_layers]))])
            general_params_dict["age.enabled"] = "true"
            general_params_dict["age.initial_value"] = start_date
            general_params_dict["isochrones.deposition_times"] = all_layers
            general_params_dict["isochrones.max_n_layers"] = 500

        outfile = f"{domain}_g{grid_resolution}m_{experiment}.nc"

        general_params_dict["output.file"] = join(dirs["state"], outfile)

        if initialstatefile is None:
            general_params_dict["bootstrap"] = ""
            general_params_dict["i"] = pism_dataname
        else:
            general_params_dict["bootstrap"] = ""
            general_params_dict["i"] = pism_dataname
            general_params_dict["input.regrid.file"] = initialstatefile
            general_params_dict["input.regrid.vars"] = regridvars

        if osize != "custom":
            general_params_dict["output.size"] = osize
        else:
            general_params_dict[
                "output.sizes.medium"
            ] = "sftgif,velsurf_mag,mask,usurf,bmelt"

        grid_params_dict = generate_grid_description(grid, domain, paleo=True)

        tlftw = 0.1

        sb_params_dict: Dict[str, Union[str, int, float]] = {
            "stress_balance.sia.bed_smoother.range": grid,
            "stress_balance.sia.enhancement_factor": combination["sia_e"],
            "stress_balance.sia.Glen_exponent": combination["sia_n"],
            "stress_balance.ssa.enhancement_factor": ssa_e,
            "stress_balance.ssa.Glen_exponent": combination["ssa_n"],
            "basal_resistance.pseudo_plastic.q": combination["pseudo_plastic_q"],
            "basal_yield_stress.mohr_coulomb.topg_to_phi.enabled": "yes",
            "basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden": combination[
                "till_effective_fraction_overburden"
            ],
            "stress_balance.blatter.enhancement_factor": combination["sia_e"],
        }
        phi_min = combination["phi_min"]
        phi_max = combination["phi_max"]
        z_min = combination["z_min"]
        z_max = combination["z_max"]

        sb_params_dict[
            "basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max"
        ] = phi_max
        sb_params_dict[
            "basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min"
        ] = phi_min
        sb_params_dict[
            "basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max"
        ] = z_max
        sb_params_dict[
            "basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min"
        ] = z_min

        sliding_law = "pseudo_plastic"
        if hasattr(combination, "sliding_law"):
            sliding_law = combination["sliding_law"]
        sb_params_dict[f"basal_resistance.{sliding_law}.enabled"] = "yes"
        stress_balance_params_dict = generate_stress_balance(
            stress_balance, sb_params_dict
        )

        pr_paleo_file_p = False
        tas_paleo_file_p = False
        atmosphere_given_file_p = False
        pr_paleo_file_p = (
            f"""$data_dir/climate/{combination["pr_paleo_file"]}"""
        )
        tas_paleo_file_p = (
            f"""$data_dir/climate/{combination["tas_paleo_file"]}"""
        )
        atmosphere_given_file_p = "$data_dir/climate/pism_SeaRISE_SMB_ext_4500m.nc"
        rho_ice = 910.0

        climate_parameters = {
            "atmosphere.given.file": atmosphere_given_file_p,
            "atmosphere.searise_greenland.file": atmosphere_given_file_p,
            "atmosphere.frac_P.file": pr_paleo_file_p,
            "atmosphere.delta_T.file": tas_paleo_file_p,
            "atmosphere.precip_scaling.file": pr_paleo_file_p,
            "surface.pdd.factor_ice": combination["f_ice"] / rho_ice,
            "surface.pdd.factor_snow": combination["f_snow"] / rho_ice,
            "surface.pdd.std_dev.value": combination["std_dev"],
            "surface.pdd.refreeze": combination["refreeze"]
        }

       
        if "atmosphere.elevation_change.temperature_lapse_rate" in combination:
            climate_parameters["atmosphere.elevation_change.temperature_lapse_rate"] = combination["atmosphere.elevation_change.temperature_lapse_rate"]
            climate_parameters["atmosphere.elevation_change.file"] = pism_dataname
        climate_params_dict = generate_climate(combination["climate"], **climate_parameters)

        hydrology_parameters = {}
        hydro_params_dict = generate_hydrology(hydrology, **hydrology_parameters)

        ocean_delta_SL_file_p = "$data_dir/ocean/pism_dSL.nc"
        ocean_parameters = {
            "sea_level.models": "constant,delta_sl",
            "ocean.delta_sl.file": ocean_delta_SL_file_p,
        }

        ocean_params_dict = generate_ocean("const", **ocean_parameters)

        calving_parameters: Dict[str, Union[str, int, float]] = {
            "calving.float_kill.calve_near_grounding_line": float_kill_calve_near_grounding_line,
            "calving.vonmises_calving.use_custom_flow_law": True,
            "calving.vonmises_calving.Glen_exponent": 3.0,
            "geometry.front_retreat.use_cfl": True,
        }
        vcm = combination["vcm"]

        try:
            vcm = float(vcm)
            calving_parameters["calving.vonmises_calving.sigma_max"] = vcm * 1e6
        except:
            vonmises_calving_threshold_file_p = "$data_dir/calving/{vcm}"
            calving_parameters[
                "calving.vonmises_calving.threshold_file"
            ] = vonmises_calving_threshold_file_p
        if "calving.eigen_calving.K" in combination:
            calving_parameters["calving.eigen_calving.K"] = combination[
                "calving.eigen_calving.K"
            ]
        if "calving.thickness_calving.threshold" in combination:
            calving_parameters["calving.thickness_calving.threshold"] = combination[
                "calving.thickness_calving.threshold"
            ]
        if "calving.thickness_calving.file" in combination:
            calving_parameters[
                "calving.thickness_calving.file"
            ] = f"""$data_dir/calving/{combination[
                "calving.thickness_calving.file"]}"""
            if "calving.thickness_calving.threshold" in calving_parameters:
                del calving_parameters["calving.thickness_calving.threshold"]

        if "calving.rate_scaling.file" in combination:
            calving_parameters[
                "calving.rate_scaling.file"
            ] = f"""$data_dir/calving/{combination[
                "calving.rate_scaling.file"]}"""
            calving_parameters["calving.rate_scaling.period"] = 0

        calving = options.calving
        calving_params_dict = generate_calving(
            calving, **calving_parameters
        )
        
        scalar_ts_dict = generate_scalar_ts(
            outfile, tsstep, odir=dirs["scalar"])
        snap_shot_dict = generate_snap_shots(outfile, "-100000, -75000, -50000, -25000,-20000,-15000,-10000,-5000,0", odir=dirs["snap"])

        solver_dict = {}


        all_params_dict = merge_dicts(
            general_params_dict,
            grid_params_dict,
            stress_balance_params_dict,
            climate_params_dict,
            ocean_params_dict,
            hydro_params_dict,
            calving_params_dict,
            scalar_ts_dict,
            snap_shot_dict,
            solver_dict,
        )

        if not spatial_ts == "none":
            exvars = spatial_ts_vars[spatial_ts]
            spatial_ts_dict = generate_spatial_ts(
                outfile, exvars, exstep, odir=dirs["spatial"]
            )

            all_params_dict = merge_dicts(all_params_dict, spatial_ts_dict)


        print("\nChecking parameters")
        print("------------------------------------------------------------")
        with xr.open_dataset(master_config_file) as ds:
            for key in all_params_dict:
                if hasattr(ds["pism_config"], key) is False:
                    print(f"  - {key} not found in pism_config")
        print("------------------------------------------------------------\n")

            
        all_params = " \\\n  ".join(
            ["-{} {}".format(k, v) for k, v in sorted(list(all_params_dict.items()))]
        )

        if commandline_options is not None:
            all_params = f"{all_params} \\\n  {commandline_options[1:-1]}"

        print("\nInput files:\n")

        check_files = []
        if pism_dataname:
            check_files.append(pism_dataname)
        if ocean_delta_SL_file_p:
            check_files.append(ocean_delta_SL_file_p)
        if atmosphere_given_file_p:
            check_files.append(atmosphere_given_file_p)
        if pr_paleo_file_p:
            check_files.append(pr_paleo_file_p)
        if tas_paleo_file_p:
            check_files.append(tas_paleo_file_p)

        for m_f in check_files:
            if "$input_dir" in m_f:
                m_f_abs = m_f.replace("$input_dir", options.input_dir)
            if "$data_dir" in m_f:
                m_f_abs = m_f.replace("$data_dir", options.data_dir)
                
            print(f"{m_f_abs}: {os.path.isfile(m_f_abs)}")
        print("\n")

        if system == "debug":
            redirect = " 2>&1 | tee {jobs}/job.${job_id}"
        else:
            redirect = " > {jobs}/job.${job_id} 2>&1"

        template = "{mpido} {pism} {params}" + redirect

        context = merge_dicts(batch_system, dirs, {"pism": pism, "params": all_params})
        cmd = template.format(**context)
        f.write(cmd)
        f.write("\n")
        f.write("\n")
        f.write(batch_system.get("footer", ""))

    scripts.append(script)


scripts = uniquify_list(scripts)
print("\n".join([script for script in scripts]))
print("\nwritten\n")
