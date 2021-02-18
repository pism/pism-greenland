#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

from datetime import datetime
import numpy as np
import pandas as pd
from scipy.io import loadmat


def to_decimal_year(date):
    year = date.year
    start_of_this_year = datetime(year=year, month=1, day=1)
    start_of_next_year = datetime(year=year + 1, month=1, day=1)
    year_elapsed = (date - start_of_this_year).total_seconds()
    year_duration = (start_of_next_year - start_of_this_year).total_seconds()
    fraction = year_elapsed / year_duration

    return date.year + fraction


if __name__ == "__main__":

    data = loadmat("xctd_fjord/JI_Mean_Each_Season.mat")
    df = pd.DataFrame(
        data=np.hstack(
            [
                data["days_since20090101"].reshape(-1, 1),
                data["depth_m"].reshape(-1, 1),
                data["temperature_mean_resolution_10m"].reshape(-1, 1),
            ]
        ),
        columns=["Time Elapsed [days since 2009-01-01]", "Depth [m]", "Temperature [Celsius]"],
    )
    df["Date"] = pd.to_datetime("2009-01-01") + pd.to_timedelta(
        df["Time Elapsed [days since 2009-01-01]"], unit="days"
    )
    time = [to_decimal_year(d) for d in df["Date"]]
    df["Year"] = time
    df = df.sort_values(by="Date")
    n = len(df)
    lon, lat = -50.27644, 69.16999
    df["Longitude [degrees_east]"] = np.repeat(lon, n).reshape(-1, 1)
    df["Latitude [degrees_north]"] = np.repeat(lat, n).reshape(-1, 1)

    df.to_csv("xctd_fjord/xctd_ilulissat_fjord.csv")
