#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

from datetime import datetime
import numpy as np
import pandas as pd
import pylab as plt
from scipy.io import loadmat


def toDecimalYear(date):
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year + 1, month=1, day=1)
    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed / yearDuration

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
        df[f"Time Elapsed [days since 2009-01-01]"], unit="days"
    )
    time = [toDecimalYear(d) for d in df["Date"]]
    df["Year"] = time
    df = df.sort_values(by="Date")
    n = len(df)
    lon, lat = -50.27644, 69.16999
    df["Longitude [degrees_east]"] = np.repeat(lon, n).reshape(-1, 1)
    df["Latitude [degrees_north]"] = np.repeat(lat, n).reshape(-1, 1)

    df.to_csv("xctd_fjord/xctd_ilulissat_fjord.csv")
