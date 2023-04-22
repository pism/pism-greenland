#!/usr/bin/env python
# Copyright (C) 2019-23 Andy Aschwanden

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
import os
import re
from glob import glob
import pandas as pd
import cftime
import xarray as xr
import time
from joblib import Parallel, delayed

import contextlib
import joblib
from tqdm import tqdm


test_file = """netcdf test_small {
  dimensions:
    bnds = 2 ;
    time = UNLIMITED ; // (12 currently)

  variables:
    double ice_mass(time) ;
      ice_mass:standard_name = "land_ice_mass" ;
      ice_mass:long_name = "mass of the ice, including seasonal cover" ;
      ice_mass:units = "kg" ;
      ice_mass:cell_methods = "time: mean" ;
      ice_mass:ancillary_variables = "ice_mass_aux" ;

    double limnsw(time) ;
      limnsw:standard_name = "land_ice_mass_not_displacing_sea_water" ;
      limnsw:long_name = "mass of the ice not displacing sea water" ;
      limnsw:units = "kg" ;
      limnsw:cell_methods = "time: mean" ;
      limnsw:ancillary_variables = "limnsw_aux" ;

    double time(time) ;
      time:standard_name = "time" ;
      time:long_name = "time" ;
      time:bounds = "time_bnds" ;
      time:units = "seconds since 2008-01-01" ;
      time:calendar = "standard" ;
      time:axis = "T" ;

    double time_bnds(time,bnds) ;

  // global attributes:
    :CDI = "Climate Data Interface version 2.1.1 (https://mpimet.mpg.de/cdi)" ;
    :Conventions = "CF-1.6" ;
    :institution = "University of Alaska Fairbanks" ;
    :proj = "epsg:3413" ;
    :frequency = "mon" ;
    :CDO = "Climate Data Operators version 2.1.1 (https://mpimet.mpg.de/cdo)" ;


  data:
    ice_mass = 2.68500643746433e+18, 2.68493466694531e+18, 2.68490395691592e+18, 2.68486729772707e+18, 2.68482531189543e+18, 2.68476719634436e+18, 2.68465977818969e+18, 2.68454523895323e+18, 2.68448060154681e+18, 2.6844324664108e+18, 2.68440161543146e+18, 2.68437852696913e+18 ;

    limnsw = 2.63391976554017e+18, 2.63394265080729e+18, 2.63395995271496e+18, 2.63397555090209e+18, 2.63398225747442e+18, 2.6339714668476e+18, 2.63391745443771e+18, 2.63385498558599e+18, 2.63383767097943e+18, 2.63383951553143e+18, 2.63385538043974e+18, 2.63387511809372e+18 ;

    time = -882273600, -879724800, -877132800, -874497600, -871862400, -869227200, -866592000, -863913600, -861278400, -858643200, -856008000, -853372800 ;

    time_bnds = 
    -883526400, -881020800, 
    -880934400, -878515200, 
    -878428800, -875836800, 
    -875750400, -873244800, 
    -873158400, -870566400, 
    -870480000, -867974400, 
    -867888000, -865296000, 
    -865209600, -862617600, 
    -862531200, -860025600, 
    -859939200, -857347200, 
    -857260800, -854755200, 
    -854668800, -852076800 ;

} // group /
"""


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


def ncfile2dataframe(infile):
    if os.path.isfile(infile):
        with xr.open_dataset(infile) as ds:
            m_id = int(re.search("id_(.+?)_", infile).group(1))
            m_dx = int(re.search("gris_g(.+?)m", infile).group(1))
            datetimeindex = ds.indexes["time"]
            nt = len(datetimeindex)
            id_S = pd.Series(data=np.repeat(m_id, nt), index=datetimeindex, name="id")
            S = [id_S]
            for m_var in ds.data_vars:
                if m_var not in ("time_bounds", "time_bnds", "timestamp"):
                    if hasattr(ds[m_var], "units"):
                        m_units = ds[m_var].units
                    else:
                        m_units = ""
                    data = np.squeeze(ds[m_var].values)
                    m_S_name = f"{m_var} ({m_units})"
                    m_S = pd.Series(data=data, index=datetimeindex, name=m_S_name)
                    S.append(m_S)
        return pd.concat(S, axis=1).reset_index()


def convert_ncfiles_to_dataframe(
    infiles, outfile: str = "scalars.parquet", format: str = "parquet", n_jobs: int = 4
):
    n_files = len(infiles)

    start_time = time.perf_counter()
    with tqdm_joblib(tqdm(desc="Processing files", total=n_files)) as progress_bar:
        result = Parallel(n_jobs=n_jobs)(
            delayed(ncfile2dataframe)(infile) for infile in infiles
        )
    finish_time = time.perf_counter()
    time_elapsed = finish_time - start_time
    print(f"Program finished in {time_elapsed:.0f} seconds")

    df = pd.concat(result)
    df.to_csv(outfile, index=False, compression="infer")
    if out_format == "csv":
        df.to_csv(outfile)
    elif out_format == "parquet":
        df.to_parquet(outfile)
    else:
        raise NotImplementedError(f"{out_format} not implemented")


if __name__ == "__main__":
    __spec__ = None

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.description = "Convert scalar netCDF variables to CSV or Parquet."
    parser.add_argument(
        "-o",
        "--output_file",
        dest="outfile",
        help="output filename",
        default="test.parquet",
    )
    parser.add_argument("-n", "--n_jobs", type=int, default=4)
    parser.add_argument("--out_format", choices=["csv", "parquet"], default="parquet")
    parser.add_argument("INFILES", nargs="*")

    options = parser.parse_args()

    outfile = options.outfile
    infiles = sorted(options.INFILES)
    out_format = options.out_format
    n_jobs = options.n_jobs

    convert_ncfiles_to_dataframe(
        infiles, outfile=outfile, format=out_format, n_jobs=n_jobs
    )
