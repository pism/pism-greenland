#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

from cftime import utime
from dateutil import rrule
from dateutil.parser import parse
from datetime import datetime
import numpy as np
import pandas as pd
import pylab as plt
from itertools import product
import os
import time

def toDecimalYear(date):
    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year + 1, month=1, day=1)
    yearElapsed = (date - startOfThisYear).total_seconds()
    yearDuration = (startOfNextYear - startOfThisYear).total_seconds()
    fraction = yearElapsed / yearDuration

    return date.year + fraction


def melting_point_temperature(depth, salinity):
    a = [-0.0575, 0.0901, -7.61e-4]
    return a[0] * salinity + a[1] + a[2] * depth


if __name__ == "__main__":

    # depths to average over
    depth_min = 225
    depth_max = 275
    # depth for freezing point calculation
    depth = 250
    salinity = 34
    grid_spacing = 4500

    calendar = "standard"
    units = "days since 1980-1-1"
    cdftime_days = utime(units, calendar)

    start_date = datetime(1980, 1, 1)
    end_date = datetime(2021, 1, 1)
    end_date_yearly = datetime(2021, 1, 2)

    # create list with dates from start_date until end_date with
    # periodicity prule.
    bnds_datelist = list(rrule.rrule(rrule.MONTHLY, dtstart=start_date, until=end_date_yearly))
    bnds_datelist_yearly = list(rrule.rrule(rrule.YEARLY, dtstart=start_date, until=end_date_yearly))

    # calculate the days since refdate, including refdate, with time being the
    bnds_interval_since_refdate = cdftime_days.date2num(bnds_datelist)
    bnds_interval_since_refdate_yearly = cdftime_days.date2num(bnds_datelist_yearly)
    time_interval_since_refdate = bnds_interval_since_refdate[0:-1] + np.diff(bnds_interval_since_refdate) / 2

    time_dict = {
        "calendar": calendar,
        "units": units,
        "time": time_interval_since_refdate,
        "time_bnds": bnds_interval_since_refdate,
    }

    dpy = np.diff(bnds_interval_since_refdate_yearly)
    dpy = np.repeat(dpy, 12, axis=0)

    decimal_time = bnds_interval_since_refdate[0:-1] / dpy + 1980

    mo_df = pd.read_csv("disko_bay_motyka.csv")

    m18_start_date_str = "2018-06-12T21:29:30"
    m18_start_date =  datetime.strptime(m18_start_date_str, '%Y-%m-%dT%H:%M:%S')
    m18_df = pd.read_csv("SBE37SMP-ODO-RS232_03715761_2018_06_12.cnv", header=454, delim_whitespace=True, usecols=[1, 2, 5], names=["Temperature", "Salinity", "Julian days"], encoding="latin1")

    m19_start_date_str = "2019-10-06T19:02:27"
    m19_start_date = datetime.strptime(m19_start_date_str, '%Y-%m-%dT%H:%M:%S')
    m19_df = pd.read_csv("sbe37smp-odo-rs232_03715761_2019_10_06.cnv", header=397, delim_whitespace=True, usecols=[0, 1, 2, 6], names=["Depth", "Temperature", "Salinity", "Time elapsed"], encoding="latin1")
    m19_time = m19_start_date +  pd.Series([pd.Timedelta(seconds=delta) for delta in m19_df["Time elapsed"]])
    m19_df["Date"] = [toDecimalYear(d) for d in m19_time]


    
    omg_df = pd.read_csv("disko_bay_omg_axctd.csv", na_values=-99.0).dropna()
    omg_df = omg_df[(omg_df["Depth"] <= depth_max) & (omg_df["Depth"] >= depth_min)]
    omg_time = pd.to_datetime(omg_df.Date, format="%m/%d/%Y %H:%M:%S")
    omg_df["Date"] = omg_time
    omg_df = omg_df.set_index("Date").drop(columns=["Unnamed: 0"])
    omg_df = omg_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    omg_time = [toDecimalYear(d) for d in omg_df.index]
    omg_df["Date"] = omg_time
    omg_df.to_csv("disko_bay_omg_depth_averaged.csv")

    ices_df = pd.read_csv("disko_bay_ices.csv")
    ices_df = ices_df[(ices_df.Depth >= depth_min) & (ices_df.Depth <= depth_max)].reset_index(drop=True)
    ices_time = pd.to_datetime(ices_df.Date, format="%Y/%m/%d %H:%M:%S")
    ices_df["Date"] = ices_time
    ices_df = ices_df.set_index("Date")
    ices_df = ices_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    ices_time = [toDecimalYear(d) for d in ices_df.index]
    ices_df["Date"] = ices_time
    ices_df.to_csv("disko_bay_ices_depth_averaged.csv")

    holl_df = pd.read_csv("disko_bay_xctd_holland.csv", parse_dates=[0])
    holl_df = holl_df[(holl_df.Depth >= depth_min) & (holl_df.Depth <= depth_max)].reset_index(drop=True)
    holl_df = holl_df.set_index("Date")
    holl_df = holl_df.groupby(pd.Grouper(freq="1D")).mean().dropna()
    holl_time = [toDecimalYear(d) for d in holl_df.index]
    holl_df["Date"] = holl_time
    holl_df.to_csv("disko_bay_ices_depth_averaged.csv")

    all_df = pd.concat([mo_df, ices_df, omg_df, holl_df])
    all_df = all_df.sort_values(by="Date")

    X_mo = mo_df.Date.values.reshape(-1, 1)
    X_ices = ices_df.Date.values.reshape(-1, 1)
    X_omg = omg_df.Date.values.reshape(-1, 1)
    X_holl = holl_df.Date.values.reshape(-1, 1)

    y_mo = mo_df.Temperature.values
    y_ices = ices_df.Temperature.values
    y_omg = omg_df.Temperature.values
    y_holl = holl_df.Temperature.values

    x_var = "Date"
    y_var = "Temperature"

    fig = plt.figure()
    ax = fig.gca()

    ax.plot(X_holl, y_holl, "o", color="#f4a582", ms=4, label="Observed (Holland)")
    ax.plot(X_mo, y_mo, "o", color="#b2182b", ms=4, label="Observed (Motyka)")
    ax.plot(X_ices, y_ices, "o", color="#92c5de", ms=4, label="Observed (ICES)")
    ax.plot(X_omg, y_omg, "o", color="#2166ac", ms=4, label="Observed (OMG)")

    ax.set_xlabel("Time")
    ax.set_ylabel("Temperature (Celsius)")
    ax.set_xlim(1980, 2021)
    ax.set_ylim(0, 5)
    plt.legend()
    fig.savefig("disko-bay-temps-obs.pdf")
    plt.cla()
    plt.close(fig)

