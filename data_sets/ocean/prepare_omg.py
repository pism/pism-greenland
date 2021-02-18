#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

from datetime import datetime
import numpy as np
import pandas as pd
import pylab as plt
from glob import glob
from netCDF4 import Dataset as NC
import os

from functools import partial
from multiprocessing import Pool


def toDecimalYear(date):
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year + 1, month=1, day=1)
    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed / yearDuration

    return date.year + fraction


def process_file(m_file):

    print(f"Processing {m_file}")

    nc = NC(m_file)

    geospatial_lat_min = float(nc.geospatial_lat_min)
    geospatial_lat_max = float(nc.geospatial_lat_max)
    geospatial_lat_units = nc.geospatial_lat_units
    geospatial_lon_min = float(nc.geospatial_lon_min)
    geospatial_lon_max = float(nc.geospatial_lon_max)
    geospatial_lon_units = nc.geospatial_lon_units

    lat = (geospatial_lat_max + geospatial_lat_min) / 2
    lon = (geospatial_lon_max + geospatial_lon_min) / 2

    density = nc.variables["density"]
    depth = nc.variables["depth"]
    temperature = nc.variables["temperature"]
    salinity = nc.variables["salinity"]

    time = nc.variables["profile_time"]
    time_units = time.units
    time_nd = np.squeeze(time[:].data)

    ref_date = nc.variables["time"]
    ref_units = ref_date.units

    start_date = pd.to_datetime(ref_units.split(" ")[-1]) + pd.to_timedelta(ref_date[0], unit="s")

    df = pd.DataFrame(
        data=np.hstack(
            [
                density[:].reshape(-1, 1),
                depth[:].reshape(-1, 1),
                temperature[:].reshape(-1, 1),
                salinity[:].reshape(-1, 1),
            ]
        ),
        columns=[
            f"Density [{density.units}]",
            f"Depth [{depth.units}]",
            f"Temperature [Celsius]",
            f"Salinity [{salinity.units}]",
        ],
    )
    df["Date"] = start_date.tz_convert(None) + pd.to_timedelta(time_nd[:], unit=time_units)
    s_time = [toDecimalYear(d) for d in df["Date"]]
    df["Year"] = s_time
    n = len(df)
    df["Longitude [degrees_east]"] = np.repeat(lon, n).reshape(-1, 1)
    df["Latitude [degrees_north]"] = np.repeat(lat, n).reshape(-1, 1)
    nc.close()

    return df


if __name__ == "__main__":
    __spec__ = None

    """
    Process OMG data from
    https://podaac-tools.jpl.nasa.gov/drive/files/allData/omg/L2/CTD/AXCTD/
    downloaded 
    $ wget -R -nc --user user --password pw -r  -np -nd  -A "*.nc" https://podaac-tools.jpl.nasa.gov/drive/files/allData/omg/L2/CTD/AXCTD 
    where user, pw are your Earthdata login username and password.

    """

    odir = "omg"
    if not os.path.isdir(odir):
        os.makedirs(odir)

    dfs = []
    files = glob("OMG_Ocean_AXCTD_L2/OMG_Ocean_AXCTD_L2_*.nc")

    n_procs = 8
    with Pool(n_procs) as pool:

        dfs = pool.map(partial(process_file), files)
        pool.close()

    df = pd.concat(dfs).reset_index(drop=True).replace(-9999, np.nan)
    df.to_csv(f"{odir}/omg_axctd_all.csv.gz", compression="gzip")

    df = df.set_index("Date")
    df = df.groupby(pd.Grouper(freq="10s")).mean().dropna()
    df["Date"] = df.index
    df.to_csv(f"{odir}/omg_axctd_all_10s_mean.csv")

    # narrow bounding box for fjord
    lon_min = -50.3
    lon_max = -50.1
    lat_min = 69.1
    lat_max = 69.3

    df_filtered = df[
        (df["Longitude [degrees_east]"] >= lon_min)
        & (df["Longitude [degrees_east]"] <= lon_max)
        & (df["Latitude [degrees_north]"] >= lat_min)
        & (df["Latitude [degrees_north]"] <= lat_max)
    ].reset_index(drop=True)
    df_filtered.to_csv(f"{odir}/omg_axctd_ilulissat_fjord_10s_mean.csv")
