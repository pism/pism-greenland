#!/usr/bin/env python
# Copyright (C) 2020 Andy Aschwanden

import numpy as np
import pandas as pd
from scipy.io import loadmat

data = loadmat("JI_Mean_Each_Season.mat")

date = pd.to_datetime(data["days_since20090101"].flatten(), unit="D", origin=pd.Timestamp("2009-01-01")).date

df = pd.DataFrame(
    data=np.vstack([date, data["depth_m"], data["temperature_mean_resolution_10m"]]).T,
    columns=["Date", "Depth", "Temperature"],
).sort_values(by=["Date", "Depth"])

df.to_csv("disko_bay_xctd_holland.csv", index=False)
