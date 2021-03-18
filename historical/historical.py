#!/usr/bin/env python
# Copyright (C) 2019-21 Andy Aschwanden

# Historical simulations for ISMIP6

import itertools
from collections import OrderedDict
import numpy as np
import os
import sys
import shlex
from os.path import join, abspath, realpath, dirname

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


grid_choices = [18000, 9000, 6000, 4500, 3600, 3000, 2400, 1800, 1500, 1200, 900, 600, 450, 300, 150, 1000]

# set up the option parser
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Generating scripts for warming experiments."
parser.add_argument("FILE", nargs=1, help="Input file to restart from", default=None)
parser.add_argument(
    "-n", "--n_procs", dest="n", type=int, help="""number of cores/processors. default=140.""", default=140
)
parser.add_argument(
    "-w", "--wall_time", dest="walltime", help="""walltime. default: 100:00:00.""", default="100:00:00"
)
parser.add_argument(
    "-q", "--queue", dest="queue", choices=list_queues(), help="""queue. default=long.""", default="long"
)
parser.add_argument("--options", dest="commandline_options", help="""Here you can add command-line options""")
parser.add_argument(
    "-d",
    "--domain",
    dest="domain",
    choices=["gris", "gris_ext", "jib", "jakobshavn", "nw", "ismip6", "qaamerujup"],
    help="sets the modeling domain",
    default="ismip6",
)
parser.add_argument("--exstep", dest="exstep", help="Writing interval for spatial time series", default="monthly")
parser.add_argument("--tsstep", dest="tsstep", help="Writing interval for scalar time series", default="daily")
parser.add_argument(
    "-f",
    "--o_format",
    dest="oformat",
    choices=["netcdf3", "netcdf4_parallel", "netcdf4_serial", "pnetcdf"],
    help="output format",
    default="netcdf4_serial",
)
parser.add_argument(
    "-L",
    "--comp_level",
    dest="compression_level",
    help="Compression level for output file. Only works with netcdf4_serial.",
    default=2,
)
parser.add_argument(
    "-g", "--grid", dest="grid", type=int, choices=grid_choices, help="horizontal grid resolution", default=9000
)
parser.add_argument("--i_dir", dest="input_dir", help="input directory", default=abspath(join(script_directory, "..")))
parser.add_argument("--o_dir", dest="output_dir", help="output directory", default="test_dir")
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
    "-b", "--bed_type", dest="bed_type", choices=list_bed_types(), help="output size type", default="wc"
)
parser.add_argument(
    "--spatial_ts",
    dest="spatial_ts",
    choices=["basic", "standard", "none", "ismip6", "strain"],
    help="output size type",
    default="standard",
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
    choices=["vonmises_calving", "hayhurst_calving"],
    help="Choose calving law",
    default="vonmises_calving",
)
parser.add_argument(
    "--stable_gl",
    dest="float_kill_calve_near_grounding_line",
    action="store_false",
    help="Stable grounding line",
    default=True,
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
    choices=["2", "3", "3a", "4", "1980", "1980a", "1980v3", "1_RAGIS"],
    help="input data set version",
    default="1_RAGIS",
)
parser.add_argument(
    "--vertical_velocity_approximation",
    dest="vertical_velocity_approximation",
    choices=["centered", "upstream"],
    help="How to approximate vertical velocities",
    default="upstream",
)
parser.add_argument("--start", help="Simulation start year", default="2008-1-1")
parser.add_argument("--end", help="Simulation end year", default="2015-1-1")
parser.add_argument(
    "-e",
    "--ensemble_file",
    dest="ensemble_file",
    help="File that has all combinations for ensemble study",
    default=None,
)

options = parser.parse_args()
commandline_options = options.commandline_options

start_date = options.start
end_date = options.end

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

spatial_ts = options.spatial_ts

bed_type = options.bed_type
calving = options.calving
exstep = options.exstep
tsstep = options.tsstep
float_kill_calve_near_grounding_line = options.float_kill_calve_near_grounding_line
grid = options.grid
hydrology = options.hydrology
stress_balance = options.stress_balance
vertical_velocity_approximation = options.vertical_velocity_approximation
version = options.version

ensemble_file = options.ensemble_file

domain = options.domain
pism_exec = generate_domain(domain)

if options.FILE is None:
    print("Missing input file")
    import sys

    sys.exit()
else:
    input_file = options.FILE[0]

if domain.lower() in ("greenland_ext", "gris_ext"):
    pism_dataname = "$input_dir/data_sets/bed_dem/pism_Greenland_ext_{}m_mcb_jpl_v{}_{}.nc".format(
        grid, version, bed_type
    )
if domain.lower() in ("ismip6"):
    pism_dataname = "$input_dir/data_sets/bed_dem/pism_Greenland_ismip6_{}m_mcb_jpl_v{}_{}.nc".format(
        grid, version, bed_type
    )
else:
    pism_dataname = "$input_dir/data_sets/bed_dem/pism_Greenland_{}m_mcb_jpl_v{}_{}.nc".format(grid, version, bed_type)

# Removed "thk" from regrid vars
# regridvars = "litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume"
regridvars = "litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk"

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

# use the actual path of the run scripts directory (we need it now and
# not during the simulation)
time_dir = join(output_dir, "time_forcing")
if not os.path.isdir(time_dir):
    os.makedirs(time_dir)

# generate the config file *after* creating the output directory
pism_config = "pism"
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

if system != "debug":
    cmd = "lfs setstripe -c -1 {}".format(dirs["output"])
    sub.call(shlex.split(cmd))
    cmd = "lfs setstripe -c -1 {}".format(dirs["spatial_tmp"])
    sub.call(shlex.split(cmd))

pism_timefile = join(time_dir, "timefile_{start}_{end}.nc".format(start=start_date, end=end_date))
try:
    os.remove(pism_timefile)
except OSError:
    pass

periodicity = "monthly"
cmd = ["create_timeline.py", "-a", start_date, "-e", end_date, "-p", periodicity, "-d", "2008-01-01", pism_timefile]
sub.call(cmd)

# ########################################################
# set up model initialization
# ########################################################

ssa_n = 3.25
ssa_e = 1.0
tefo = 0.020
phi_min = 5.0
phi_max = 40.0
topg_min = -700
topg_max = 700

combinations = np.genfromtxt(ensemble_file, dtype=None, encoding=None, delimiter=",", skip_header=1)

scripts = []
scripts_post = []

simulation_start_year = options.start
simulation_end_year = options.end

batch_header, batch_system = make_batch_header(system, nn, walltime, queue)
post_header = make_batch_post_header(system)

for n, combination in enumerate(combinations):

    (
        run_id,
        climate,
        hydrology,
        frontal_melt,
        climate_file,
        runoff_file,
        frontal_melt_file,
        mbp_file,
        fm_a,
        fm_b,
        fm_alpha,
        fm_beta,
        vcm,
        calving_threshold,
        gamma_T,
        salinity,
        ppq,
        sia_e,
    ) = combination

    ttphi = "{},{},{},{}".format(phi_min, phi_max, topg_min, topg_max)

    name_options = {}
    try:
        name_options["id"] = "{:03d}".format(int(run_id))
    except:
        name_options["id"] = "{}".format(run_id)

    vversion = "v" + str(version)
    full_exp_name = "_".join([vversion, "_".join(["_".join([k, str(v)]) for k, v in list(name_options.items())])])
    full_outfile = "g{grid}m_{experiment}.nc".format(grid=grid, experiment=full_exp_name)

    experiment = "_".join(
        [
            vversion,
            "_".join(["_".join([k, str(v)]) for k, v in list(name_options.items())]),
            "{}".format(start_date),
            "{}".format(end_date),
        ]
    )

    script = join(scripts_dir, "{}_g{}m_{}.sh".format(domain, grid, experiment))
    scripts.append(script)

    for filename in script:
        try:
            os.remove(filename)
        except OSError:
            pass

    with open(script, "w") as f:

        f.write(batch_header)
        f.write(run_header)

        outfile = "{domain}_g{grid}m_{experiment}.nc".format(domain=domain.lower(), grid=grid, experiment=experiment)

        pism = generate_prefix_str(pism_exec)

        general_params_dict = {
            "profile": join(dirs["performance"], "profile_${job_id}.py".format(**batch_system)),
            "time_file": pism_timefile,
            "o": join(dirs["state"], outfile),
            "o_format": oformat,
            "output.compression_level": compression_level,
            "config_override": "$config",
        }

        general_params_dict["bootstrap"] = ""
        general_params_dict["i"] = pism_dataname
        general_params_dict["regrid_file"] = input_file
        general_params_dict["regrid_vars"] = regridvars

        if osize != "custom":
            general_params_dict["o_size"] = osize
        else:
            general_params_dict["output.sizes.medium"] = "sftgif,velsurf_mag,mask,usurf"

        grid_params_dict = generate_grid_description(grid, domain)

        THRESHOLD = 4.5e4  #  stress threshold
        FRACRATE = 0.5  #  fracture rate
        HEALTHRESHOLD = 2.0e-10  #  healing threshold
        HEALRATE = 0.05  #  healing rate
        SOFTRES = 0.01  #  softening residual (avoid viscosity from degeneration), value 1 inhibits softening effect

        sb_params_dict = {
            "sia_e": sia_e,
            "ssa_e": ssa_e,
            "ssa_n": ssa_n,
            "pseudo_plastic_q": ppq,
            "till_effective_fraction_overburden": tefo,
            "vertical_velocity_approximation": vertical_velocity_approximation,
            "stress_balance.blatter.enhancement_factor": sia_e
            # "fractures": True,
            # "fracture_parameters": "{},{},{},{}".format(FRACRATE, THRESHOLD, HEALRATE, HEALTHRESHOLD),
            # "write_fd_fields": True,
            # "scheme_fd2d": True,
        }

        sb_params_dict["topg_to_phi"] = ttphi

        stress_balance_params_dict = generate_stress_balance(stress_balance, sb_params_dict)

        climate_parameters = {
            "climate_forcing.buffer_size": 367,
            "surface_given_file": "$input_dir/data_sets/ismip6/{}".format(climate_file),
        }

        climate_params_dict = generate_climate(climate, **climate_parameters)

        hydrology_parameters = {
            "hydrology.routing.include_floating_ice": True,
            "hydrology.surface_input_file": "$input_dir/data_sets/ismip6/{}".format(runoff_file),
            "hydrology.routing.add_water_input_to_till_storage": False,
            "hydrology.add_water_input_to_till_storage": False,
        }

        hydro_params_dict = generate_hydrology(hydrology, **hydrology_parameters)

        if frontal_melt == "discharge_routing":
            hydrology_parameters["hydrology.surface_input.file"] = "$input_dir/data_sets/ismip6/{}".format(
                frontal_melt_file
            )

            frontalmelt_parameters = {
                "frontal_melt": "routing",
                "frontal_melt.routing.file": "$input_dir/data_sets/ocean/{}".format(frontal_melt_file),
            }
        else:
            frontalmelt_parameters = {
                "frontal_melt": "discharge_given",
                "frontal_melt.discharge_given.file": "$input_dir/data_sets/ocean/{}".format(frontal_melt_file),
            }

        frontalmelt_params_dict = frontalmelt_parameters

        # Need to add salinity first
        ocean_parameters = {
            "ocean.th.file": "$input_dir/data_sets/ocean/{}".format(frontal_melt_file),
            "ocean.th.clip_salinity": False,
            "ocean.th.gamma_T": gamma_T,
        }
        if salinity:
            ocean_parameters["constants.sea_water.salinity"] = salinity

        if mbp_file:
            ocean_parameters["ocean.frac_MBP.file"] = "$input_dir/data_sets/melange/{}".format(mbp_file)
            ocean_parameters["ocean.frac_MBP.period"] = 1
            ocean_parameters["ocean.melange_back_pressure_fraction"] = 1.0
            ocean_params_dict = generate_ocean("th_mbp", **ocean_parameters)
        else:
            ocean_params_dict = generate_ocean("th", **ocean_parameters)

        try:
            vcm = float(vcm)
            calving_parameters = {
                "float_kill_calve_near_grounding_line": float_kill_calve_near_grounding_line,
                "calving.vonmises_calving.sigma_max": vcm * 1e6,
                "calving.vonmises_calving.use_custom_flow_law": True,
                "calving.vonmises_calving.Glen_exponent": 3.0,
                "calving.thickness_calving.threshold": calving_threshold,
            }
        except:
            calving_parameters = {
                "float_kill_calve_near_grounding_line": float_kill_calve_near_grounding_line,
                "calving.vonmises_calving.threshold_file": f"$input_dir/data_sets/calving/{vcm}"
                "calving.vonmises_calving.use_custom_flow_law": True,
                "calving.vonmises_calving.Glen_exponent": 3.0,
            }
        calving_params_dict = generate_calving(calving, **calving_parameters)

        scalar_ts_dict = generate_scalar_ts(outfile, tsstep, odir=dirs["scalar"])

        solver_dict = {}

        all_params_dict = merge_dicts(
            general_params_dict,
            grid_params_dict,
            stress_balance_params_dict,
            climate_params_dict,
            ocean_params_dict,
            hydro_params_dict,
            frontalmelt_params_dict,
            calving_params_dict,
            scalar_ts_dict,
            solver_dict,
        )

        if not spatial_ts == "none":
            exvars = spatial_ts_vars[spatial_ts]
            spatial_ts_dict = generate_spatial_ts(outfile, exvars, exstep, odir=dirs["spatial_tmp"], split=False)

            all_params_dict = merge_dicts(all_params_dict, spatial_ts_dict)

        all_params = " \\\n  ".join(["-{} {}".format(k, v) for k, v in list(all_params_dict.items())])
        if commandline_options is not None:
            all_params = f"{all_params} \\\n  {commandline_options[1:-1]}"

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
        if not spatial_ts == "none":
            f.write(
                "mv {tmpfile} {ofile}\n".format(
                    tmpfile=spatial_ts_dict["extra_file"], ofile=join(dirs["spatial"], "ex_" + outfile)
                )
            )
        f.write("\n")
        f.write(batch_system.get("footer", ""))

    scripts.append(script)


scripts = uniquify_list(scripts)
print("\n".join([script for script in scripts]))
print("\nwritten\n")
