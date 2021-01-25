#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

import numpy as np
import pandas as pd
import pylab as plt
from glob import glob

if __name__ == "__main__":

    dfs = []
    files = glob("omg_axctd/OMG_*.edf")
    for m_file in files:
        print(f"Processing {m_file}")
        header_info = {"Date": [], "Time": [], "Longitude": [], "Latitude": []}
        with open(m_file, encoding="latin-1") as f:
            for k, line in enumerate(f.readlines(1500)):
                for it in header_info.keys():
                    if line.find(it) == 0:
                        header_info[it].append(k)
        with open(m_file, encoding="latin-1") as f:
            data = f.readlines(1500)
            date = data[header_info["Date"][0]].split(" : ")[-1].replace("\n", " ").strip()
            time = data[header_info["Time"][0]].split(" : ")[-1].replace("\n", " ").strip()
            if "N\n" in data[header_info["Latitude"][0]]:
                degrees, minutes = (
                    data[header_info["Latitude"][0]].split(":")[-1].lstrip().replace("N\n", "").split(" ")
                )
                lat = float(degrees) + float(minutes) / 60
            elif "S\n" in data[header_info["Latitude"][0]]:
                degrees, minutes = (
                    data[header_info["Latitude"][0]].split(":")[-1].lstrip().replace("S\n", "").split(" ")
                )
                lat = -(float(degrees) + float(minutes) / 60)
            else:
                lat = data[header_info["Latitude"][0]].split(":")[-1].replace("\n", " ").strip()
            if "W\n" in data[header_info["Longitude"][0]]:
                degrees, minutes = (
                    data[header_info["Longitude"][0]].split(":")[-1].lstrip().replace("W\n", "").split(" ")
                )
                lon = -(float(degrees) + float(minutes) / 60)
            elif "E\n" in data[header_info["Longitude"][0]]:
                degrees, minutes = (
                    data[header_info["Longitude"][0]].split(":")[-1].lstrip().replace("E\n", "").split(" ")
                )
                lon = float(degrees) + float(minutes) / 60
            else:
                lon = data[header_info["Longitude"][0]].split(":")[-1].replace("\n", " ").strip()

        s_df = pd.read_csv(
            m_file,
            sep="\t",
            header=49,
            usecols=[2, 3, 5],
            encoding="latin-1",
            names=["depth", "temperature", "salinity"],
        )
        n = len(s_df.depth)
        if n != 0:
            m_df = pd.DataFrame(
                data=np.hstack(
                    [
                        np.repeat(date + " " + time, n).reshape(-1, 1),
                        np.repeat(lon, n).reshape(-1, 1),
                        np.repeat(lat, n).reshape(-1, 1),
                        s_df.depth.values[0:n].reshape(-1, 1),
                        s_df.temperature.values[0:n].reshape(-1, 1),
                        s_df.salinity.values[0:n].reshape(-1, 1),
                    ]
                ),
                columns=["Date", "Longitude", "Latitude", "Depth", "Temperature", "Salinity"],
            )
            dfs.append(m_df)

    df = pd.concat(dfs).reset_index(drop=True)
    df = df.astype({"Longitude": float, "Latitude": float, "Depth": float, "Temperature": float})
    df.to_csv("omg_axctd.csv")

    # depths to average over
    depth_min = 225
    depth_max = 275

    # bounding box by Khazendar (2019)
    lon_min = -54.0208
    lon_max = -52.0096
    lat_min = 68.8608
    lat_max = 69.321

    df_da = df[(df.Depth >= depth_min) & (df.Depth <= depth_max)].reset_index(drop=True)
    df_da.to_csv("omg_axctd_da_225_275.csv")

    df_filtered = df[
        (df.Longitude >= lon_min) & (df.Longitude <= lon_max) & (df.Latitude >= lat_min) & (df.Latitude <= lat_max)
    ].reset_index(drop=True)
    df_filtered.to_csv("omg_axctd_disko_bay_khazendar.csv")

    # narrow bounding box for fjord
    lon_min = -50.3
    lon_max = -50.1
    lat_min = 69.1
    lat_max = 69.3

    df_filtered = df[
        (df.Longitude >= lon_min) & (df.Longitude <= lon_max) & (df.Latitude >= lat_min) & (df.Latitude <= lat_max)
    ].reset_index(drop=True)
    df_filtered.to_csv(
        "disk_bay_fjord_omg_axctd.csv", columns=["Date", "Longitude", "Latitude", "Depth", "Temperature", "Salinity"]
    )

    # narrow bounding box for fjord
    lon_min = -54.0208
    lon_max = -50.1
    lat_min = 68.8608
    lat_max = 69.321

    df_filtered = df[
        (df.Longitude >= lon_min) & (df.Longitude <= lon_max) & (df.Latitude >= lat_min) & (df.Latitude <= lat_max)
    ].reset_index(drop=True)
    df_filtered.to_csv(
        "disko_bay_omg_axctd.csv", columns=["Date", "Longitude", "Latitude", "Depth", "Temperature", "Salinity"]
    )
