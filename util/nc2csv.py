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


def process_file(infile):
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
    return pd.concat(S, axis=1).reset_index()


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
    idx_norm_years = 0
    n_files = len(infiles)

    start_time = time.perf_counter()
    with tqdm_joblib(tqdm(desc="Processing file", total=n_files)) as progress_bar:
        result = Parallel(n_jobs=n_jobs)(
            delayed(process_file)(infile) for infile in infiles
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
