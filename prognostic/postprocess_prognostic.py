#!/usr/bin/env python
# Copyright (C) 2019 Andy Aschwanden

import os
from os.path import abspath, basename, join
import glob
import numpy as np
import re
from dateutil import rrule
from dateutil.parser import parse

from argparse import ArgumentParser

import multiprocessing

# multiprocessing.set_start_method("forkserver", force=True)
from multiprocessing import Pool
from functools import partial

from netCDF4 import Dataset as NC
from cftime import utime, datetime
from cdo import Cdo

cdo = Cdo()


def adjust_timeline(
    filename,
    start_date="2015-1-1",
    interval=1,
    interval_type="mid",
    bounds=True,
    periodicity="yearly".upper(),
    ref_date="2008-1-1",
    ref_unit="days",
    calendar="standard",
):
    nc = NC(filename, "a")
    nt = len(nc.variables["time"])

    time_units = "%s since %s" % (ref_unit, ref_date)
    time_calendar = calendar

    cdftime = utime(time_units, time_calendar)

    # create a dictionary so that we can supply the periodicity as a
    # command-line argument.
    pdict = {}
    pdict["SECONDLY"] = rrule.SECONDLY
    pdict["MINUTELY"] = rrule.MINUTELY
    pdict["HOURLY"] = rrule.HOURLY
    pdict["DAILY"] = rrule.DAILY
    pdict["WEEKLY"] = rrule.WEEKLY
    pdict["MONTHLY"] = rrule.MONTHLY
    pdict["YEARLY"] = rrule.YEARLY
    prule = pdict[periodicity]

    # reference date from command-line argument
    r = time_units.split(" ")[2].split("-")
    refdate = datetime(int(r[0]), int(r[1]), int(r[2]))

    # create list with dates from start_date for nt counts
    # periodicity prule.
    bnds_datelist = list(rrule.rrule(freq=prule, dtstart=parse(start_date), count=nt + 1, interval=interval))

    # calculate the days since refdate, including refdate, with time being the
    bnds_interval_since_refdate = cdftime.date2num(bnds_datelist)
    if interval_type == "mid":
        # mid-point value:
        # time[n] = (bnds[n] + bnds[n+1]) / 2
        time_interval_since_refdate = bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2
    elif interval_type == "start":
        time_interval_since_refdate = bnds_interval_since_refdate[:-1]
    else:
        time_interval_since_refdate = bnds_interval_since_refdate[1:]

    # create a new dimension for bounds only if it does not yet exist
    time_dim = "time"
    if time_dim not in list(nc.dimensions.keys()):
        nc.createDimension(time_dim)

    # variable names consistent with PISM
    time_var_name = "time"
    bnds_var_name = "time_bnds"

    # create time variable
    if time_var_name not in nc.variables:
        time_var = nc.createVariable(time_var_name, "d", dimensions=(time_dim))
    else:
        time_var = nc.variables[time_var_name]
    time_var[:] = time_interval_since_refdate
    time_var.units = time_units
    time_var.calendar = time_calendar
    time_var.standard_name = time_var_name
    time_var.axis = "T"

    if bounds:
        # create a new dimension for bounds only if it does not yet exist
        bnds_dim = "nb2"
        if bnds_dim not in list(nc.dimensions.keys()):
            nc.createDimension(bnds_dim, 2)

        # create time bounds variable
        if bnds_var_name not in nc.variables:
            time_bnds_var = nc.createVariable(bnds_var_name, "d", dimensions=(time_dim, bnds_dim))
        else:
            time_bnds_var = nc.variables[bnds_var_name]
        time_bnds_var[:, 0] = bnds_interval_since_refdate[0:-1]
        time_bnds_var[:, 1] = bnds_interval_since_refdate[1::]
        time_var.bounds = bnds_var_name
    else:
        delattr(time_var, "bounds")

    nc.close()


def get_var(m_string):
    """
    Return the variable name
    """
    m_var = None
    for ismip6var in ISMIP6.keys():
        if re.search(ismip6var, m_string) is not None:
            m_var = ismip6var

    return m_var


def process_file(a_file, metadata):

    m_file = basename(a_file)
    m_exp = m_file.split("_")[-1].split(".")[0]
    m_var = get_var(m_file)

    EXP_GRID = "_".join([m_exp, "01"])

    base_dir = metadata["base_dir"]

    if m_var is not None:

        print("Processing {}".format(m_file))
        project_dir = os.path.join(base_dir, GROUP, MODEL, EXP_GRID)
        if not os.path.exists(project_dir):
            os.makedirs(project_dir)

        o_file = join(project_dir, m_file)

        if ISMIP6[m_var]["dims"] == str(1):
            time_interval = 1
        elif ISMIP6[m_var]["dims"] == str(2):
            time_interval = 5
        else:
            print("Wrong dims for time interval")

        if ISMIP6[m_var]["type"] == "state":
            print("  Saving {}".format(o_file))
            cdo.seltimestep(
                "1/1000", input="-setmissval,1.e20 {}".format(a_file), output=o_file, options="-f nc4 -z zip_3 -O -L"
            )
            adjust_timeline(o_file, interval=time_interval, interval_type="start", bounds=False)
        elif ISMIP6[m_var]["type"] == "flux":
            print("  Saving {}".format(o_file))
            cdo.seltimestep(
                "2/1000", input="-setmissval,1.e20 {}".format(a_file), output=o_file, options="-f nc4 -z zip_3 -O -L"
            )
            adjust_timeline(o_file, interval=time_interval, interval_type="mid", bounds=True)
        else:
            print("how did I get here")


# from ISMIP6 import *

# Set up the option parser
parser = ArgumentParser()
parser.description = "Script to make ISMIP6-conforming time series."
parser.add_argument("INDIR", nargs=1)
parser.add_argument("--model", dest="model", type=str, help="""Model ID""", default="1")
parser.add_argument("-o", dest="base_dir", type=str, help="""Basedirectory for output""", default=".")
parser.add_argument(
    "--resource_dir", dest="resource_dir", type=str, help="""Directory with ISMIP6 resources""", default="."
)
parser.add_argument(
    "-n", "--n_procs", dest="n_procs", type=int, help="""number of cores/processors. default=4.""", default=4
)

options = parser.parse_args()
n_procs = options.n_procs
in_dir = abspath(options.INDIR[0])
base_dir = abspath(options.base_dir)

scalar_dir = join(in_dir, "scalar_processed")
spatial_dir = join(in_dir, "spatial")
IS = "GIS"
GROUP = "UAF"
MODEL = options.model
project = "{IS}_{GROUP}_{MODEL}".format(IS=IS, GROUP=GROUP, MODEL=MODEL)

try:
    ismip6resources = np.genfromtxt(
        "../resources/ismip6vars.csv", dtype=None, encoding=None, delimiter=",", autostrip=True
    )
except:
    ismip6resources = np.genfromtxt("../resources/ismip6vars.csv", dtype=None, delimiter=",", autostrip=True)

ISMIP6 = {}
keys = ismip6resources[0]
for v in ismip6resources[1::]:
    ISMIP6[v[0]] = dict(zip(keys[1::], v[1::]))


if __name__ == "__main__":

    print("Postprocessing PISM simulation for ISMIP6")
    print("-------------------------------------------------\n")
    print("indir: {}".format(in_dir))
    print("basedir: {}".format(base_dir))

    files = []
    files.extend(glob.glob(join(scalar_dir, "*.nc")))
    files.extend(glob.glob(join(spatial_dir, "*.nc")))

    metadata = {"base_dir": base_dir}
    pool = Pool(n_procs)
    pool.map(partial(process_file, metadata=metadata), files)
    pool.terminate()
