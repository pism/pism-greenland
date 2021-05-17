#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

import pandas as pd
from datetime import datetime


def to_decimal_year(date):
    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction


# depths to average over
depth_min = 225
depth_max = 275

if __name__ == "__main__":

    odir = "ices"
    df_b = pd.read_csv("ices/ices_data_bottle.csv.gz", low_memory=False)
    df_ctd = pd.read_csv("ices/ices_data_ctd.csv.gz", low_memory=False)
    df = pd.concat([df_b, df_ctd]).reset_index(drop=True)
    df = df.rename(
        columns={
            "TEMP [deg C]": "Temperature [Celsius]",
            "PRES [db]": "Depth [m]",
            "Longitude [degrees_east]": "Longitude [degrees_east]",
            "Latitude [degrees_north]": "Latitude [degrees_north]",
            "PSAL [psu]": "Salinity [g/kg]",
        }
    ).reset_index(drop=True)
    df["Date"] = pd.to_datetime(df["yyyy-mm-ddThh:mm"])

    df = df.astype(
        {
            "Longitude [degrees_east]": float,
            "Latitude [degrees_north]": float,
            "Depth [m]": float,
            "Temperature [Celsius]": float,
        }
    )
    df = df[
        [
            "Date",
            "Longitude [degrees_east]",
            "Latitude [degrees_north]",
            "Depth [m]",
            "Temperature [Celsius]",
            "Salinity [g/kg]",
        ]
    ].dropna()
    time = [to_decimal_year(d) for d in df["Date"]]
    df["Year"] = time

    df.to_csv(f"{odir}/ices_all.csv.gz", compression="gzip")

    lon_min = -52.75
    lon_max = -51.05
    lat_min = 68.50
    lat_max = 69.50

    df_filtered = df[
        (df["Longitude [degrees_east]"] >= lon_min)
        & (df["Longitude [degrees_east]"] <= lon_max)
        & (df["Latitude [degrees_north]"] >= lat_min)
        & (df["Latitude [degrees_north]"] <= lat_max)
    ].reset_index(drop=True)
    df_filtered.to_csv(f"{odir}/ices_disko_bay.csv")
    df_filtered = df_filtered[(df_filtered["Depth [m]"] <= depth_max) & (df_filtered["Depth [m]"] >= depth_min)]
    df_filtered.to_csv(f"{odir}/ices_disko_bay_250m.csv")
