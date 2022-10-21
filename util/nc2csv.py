#!/usr/bin/env python
# Copyright (C) 2019-21 Andy Aschwanden

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
from netCDF4 import Dataset as NC
import os
import re
from glob import glob
import pandas as pd
import cftime
import xarray as xr

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Convert scalar netCDF variables to CSV."
parser.add_argument(
    "-o",
    "--output_file",
    dest="outfile",
    help="Output CSV file name",
    default="test.csv.gz",
)
parser.add_argument("INFILES", nargs="*")

options = parser.parse_args()

outfile = options.outfile
infiles = options.INFILES
idx_norm_years = 0


dfs = []
for infile in infiles:
    print(f"Processing {infile}")
    if os.path.isfile(infile):
        ds = xr.open_dataset(infile)
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
    dfs.append(pd.concat(S, axis=1).reset_index())
df = pd.concat(dfs)
df.to_csv(outfile, index=False, compression="infer")
