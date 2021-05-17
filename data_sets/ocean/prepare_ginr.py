#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

from datetime import datetime, timedelta
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


def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    return datetime.fromordinal(int(datenum)) + timedelta(days=days) - timedelta(days=366)


rho_sea_water = 1.0273
# depths to average over
depth_min = 225
depth_max = 275

if __name__ == "__main__":

    m_dfs = []

    data = loadmat("ginr/GINR_2007_CTD.mat", squeeze_me=True, simplify_cells=True)
    dfs = []
    for cast in data["GINR_2007_CTD"]:
        print(f"Processing {cast}")
        c_data = data["GINR_2007_CTD"][cast]
        n = len(c_data["temp"])
        lat = c_data["lat"]
        lon = c_data["lon"]
        start_time = c_data["startTime"]

        date = pd.to_datetime(datenum_to_datetime(start_time))
        df = pd.DataFrame(
            data=np.hstack(
                [
                    np.repeat(date, n).reshape(-1, 1),
                    c_data["pres_db"].reshape(-1, 1),
                    c_data["pres_db"].reshape(-1, 1) / rho_sea_water,
                    c_data["salt_abs"].reshape(-1, 1),
                    c_data["temp"].reshape(-1, 1),
                    c_data["pot_temp"].reshape(-1, 1),
                    c_data["cond_Sm"].reshape(-1, 1),
                    np.repeat(lon, n).reshape(-1, 1),
                    np.repeat(lat, n).reshape(-1, 1),
                    np.repeat(cast, n).reshape(-1, 1),
                ]
            ),
            columns=[
                "Date",
                "Pressure [db]",
                "Depth [m]",
                "Salinity [g/kg]",
                "Temperature [Celsius]",
                "Potential Temperature [Celsius]",
                "Conductivity []",
                "Longitude [degrees_east]",
                "Latitude [degrees_north]",
                "Cast",
            ],
        )
        dfs.append(df)

    df = pd.concat(dfs)
    m_dfs.append(df)

    data = loadmat("ginr/GINR_91_94.mat", squeeze_me=True, simplify_cells=True)
    dfs = []
    for c_data in data["trawl_data"]:
        cast = c_data["metadata"]["fname"].split("/")[-1].split(".")[0]
        print(f"Processing {cast}")
        try:
            n = len(c_data["data"]["temp"])
        except:
            n = 1
        lat = c_data["lat"]
        lon = c_data["lon"]
        start_time = c_data["datenum"]
        date = pd.to_datetime(datenum_to_datetime(start_time))
        if n != 1:
            df = pd.DataFrame(
                data=np.hstack(
                    [
                        np.repeat(date, n).reshape(-1, 1),
                        c_data["data"]["prdM"].reshape(-1, 1),
                        c_data["data"]["prdM"].reshape(-1, 1) / rho_sea_water,
                        c_data["data"]["temp"].reshape(-1, 1),
                        c_data["data"]["sal00"].reshape(-1, 1),
                        np.repeat(lon, n).reshape(-1, 1),
                        np.repeat(lat, n).reshape(-1, 1),
                        np.repeat(cast, n).reshape(-1, 1),
                    ]
                ),
                columns=[
                    "Date",
                    "Pressure [db]",
                    "Depth [m]",
                    "Temperature [Celsius]",
                    "Salinity [g/kg]",
                    "Longitude [degrees_east]",
                    "Latitude [degrees_north]",
                    "Cast",
                ],
            )
        else:
            df = pd.DataFrame(
                data=np.hstack(
                    [
                        np.array([date]).reshape(-1, 1),
                        np.array([c_data["data"]["prdM"]]).reshape(-1, 1),
                        np.array([c_data["data"]["prdM"]]).reshape(-1, 1) / rho_sea_water,
                        np.array([c_data["data"]["temp"]]).reshape(-1, 1),
                        np.array([c_data["data"]["sal00"]]).reshape(-1, 1),
                        np.array([lon]).reshape(-1, 1),
                        np.array([lat]).reshape(-1, 1),
                        np.array([cast]).reshape(-1, 1),
                    ]
                ),
                columns=[
                    "Date",
                    "Pressure [db]",
                    "Depth [m]",
                    "Temperature [Celsius]",
                    "Salinity [g/kg]",
                    "Longitude [degrees_east]",
                    "Latitude [degrees_north]",
                    "Cast",
                ],
            )
        dfs.append(df)

    df = pd.concat(dfs)
    m_dfs.append(df)

    s_dfs = []
    for m_file in ["ginr/shrimp_data_2005_2010.mat", "ginr/shrimp_data_2011.mat"]:
        data = loadmat(m_file, squeeze_me=True, simplify_cells=True)
        dfs = []
        c_data = data["shrimp_data"]
        n = len(c_data["year1"])
        lat = (np.array(c_data["lat2"]) + np.array(c_data["lat1"])) / 2
        lon = (np.array(c_data["lon2"]) + np.array(c_data["lon1"])) / 2
        year = (np.array(c_data["year2"]) + np.array(c_data["year1"])) / 2
        year = year * (year != 11) + (year + 2000) * (year == 11)
        month = (np.array(c_data["month2"]) + np.array(c_data["month1"])) / 2
        day = (np.array(c_data["day2"]) + np.array(c_data["day1"])) / 2
        depth = np.array((c_data["depth2"]) + np.array(c_data["depth1"])) / 2
        temp = np.array(c_data["btm_temp"])
        date = pd.to_datetime(pd.DataFrame({"year": year, "month": month, "day": day}))
        c1 = np.array(c_data["trip_id"])
        c2 = np.array(c_data["station_id"])
        cast = np.array([f"Station {c2_i} Trip {c1_i}" for c1_i, c2_i in zip(c1, c2)])

        date = pd.to_datetime(datenum_to_datetime(start_time))
        df = pd.DataFrame(
            data=np.hstack(
                [
                    depth.reshape(-1, 1) * rho_sea_water,
                    depth.reshape(-1, 1),
                    temp.reshape(-1, 1),
                    lon.reshape(-1, 1),
                    lat.reshape(-1, 1),
                    cast.reshape(-1, 1),
                ]
            ),
            columns=[
                "Pressure [db]",
                "Depth [m]",
                "Temperature [Celsius]",
                "Longitude [degrees_east]",
                "Latitude [degrees_north]",
                "Cast",
            ],
        )
        df["Date"] = date
        s_dfs.append(df)
    df = pd.concat(s_dfs)
    m_dfs.append(df)

    df = pd.concat(m_dfs)
    time = [to_decimal_year(d) for d in df["Date"]]
    df["Year"] = time
    df = df.astype(
        {
            "Pressure [db]": float,
            "Depth [m]": float,
            "Salinity [g/kg]": float,
            "Temperature [Celsius]": float,
            "Longitude [degrees_east]": float,
            "Latitude [degrees_north]": float,
        }
    )
    df.to_csv("ginr/ginr_all.csv")

    lon_min = -52.75
    lon_max = -51.05
    lat_min = 68.50
    lat_max = 69.50

    df = df[
        (df["Longitude [degrees_east]"] >= lon_min)
        & (df["Longitude [degrees_east]"] <= lon_max)
        & (df["Latitude [degrees_north]"] >= lat_min)
        & (df["Latitude [degrees_north]"] <= lat_max)
    ].reset_index(drop=True)
    df.to_csv("ginr/ginr_disko_bay.csv")

    df = df[(df["Depth [m]"] <= depth_max) & (df["Depth [m]"] >= depth_min)]
    df.to_csv("ginr/ginr_disko_bay_250m.csv")
