#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

from datetime import datetime
import numpy as np
import pandas as pd
from glob import glob
from netCDF4 import Dataset as NC
import os

from functools import partial
from multiprocessing import Pool


def to_decimal_year(date):
    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction


def process_file(m_file):

    print(f"Processing {m_file}")

    nc = NC(m_file)

    geospatial_lat_min = float(nc.geospatial_lat_min)
    geospatial_lat_max = float(nc.geospatial_lat_max)
    geospatial_lon_min = float(nc.geospatial_lon_min)
    geospatial_lon_max = float(nc.geospatial_lon_max)

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

    start_date = pd.to_datetime(ref_units.split(" ")[-1]) + pd.to_timedelta(
        ref_date[0], unit="s"
    )

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
            "Depth [m]",
            "Temperature [Celsius]",
            "Salinity [g/kg]",
        ],
    )
    df["Date"] = start_date.tz_convert(None) + pd.to_timedelta(
        time_nd[:], unit=time_units
    )
    s_time = [to_decimal_year(d) for d in df["Date"]]
    df["Year"] = s_time
    n = len(df)
    df["Longitude [degrees_east]"] = np.repeat(lon, n).reshape(-1, 1)
    df["Latitude [degrees_north]"] = np.repeat(lat, n).reshape(-1, 1)
    nc.close()

    return df


rho_sea_water = 1.0273
# depths to average over
depth_min = 225
depth_max = 275

if __name__ == "__main__":
    __spec__ = None

    """
    Process OMG data from
    https://podaac-tools.jpl.nasa.gov/drive/files/allData/omg/L2/CTD/AXCTD/
    downloaded
    $ wget -R -nc --user user --password pw -r  -np -nd  -A "*.nc" \
      https://podaac-tools.jpl.nasa.gov/drive/files/allData/omg/L2/CTD/AXCTD
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
    time = [to_decimal_year(d) for d in df["Date"]]
    df["Year"] = time

    df.to_csv(f"{odir}/omg_axctd_all.csv.gz", compression="gzip")

    df = df.set_index("Date")
    df = df.groupby(pd.Grouper(freq="10s")).mean().dropna()
    df["Date"] = df.index
    df.to_csv(f"{odir}/omg_axctd_all_10s_mean.csv")

    lon_min = -52.75
    lon_max = -51.05
    lat_min = 68.50
    lat_max = 69.50

    df_bay = df[
        (df["Longitude [degrees_east]"] >= lon_min)
        & (df["Longitude [degrees_east]"] <= lon_max)
        & (df["Latitude [degrees_north]"] >= lat_min)
        & (df["Latitude [degrees_north]"] <= lat_max)
    ].reset_index(drop=True)
    df_bay.to_csv(f"{odir}/omg_axctd_disko_bay_10s_mean.csv")

    df_bay = df_bay[
        (df_bay["Depth [m]"] <= depth_max) & (df_bay["Depth [m]"] >= depth_min)
    ]
    df_bay.to_csv(f"{odir}/omg_axctd_disko_bay_10s_mean_250m.csv")

    # narrow bounding box for fjord
    lon_min = -51.00
    lon_max = -49.50
    lat_min = 69.00
    lat_max = 69.35

    df_fjord = df[
        (df["Longitude [degrees_east]"] >= lon_min)
        & (df["Longitude [degrees_east]"] <= lon_max)
        & (df["Latitude [degrees_north]"] >= lat_min)
        & (df["Latitude [degrees_north]"] <= lat_max)
    ].reset_index(drop=True)
    df_fjord.to_csv(f"{odir}/omg_axctd_ilulissat_fjord_10s_mean.csv")

    df_fjord = df_fjord[
        (df_fjord["Depth [m]"] <= depth_max) & (df_fjord["Depth [m]"] >= depth_min)
    ]
    df_fjord.to_csv(f"{odir}/omg_axctd_ilulissat_fjord_10s_mean_250m.csv")
