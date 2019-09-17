#!/usr/bin/env python
# Copyright (C) 2019 Andy Aschwanden

import pandas as pd
import os

csv_files = [
    "../../uncertainty_qunatification/prognostic_core.csv",
    "../../uncertainty_qunatification/prognostic_open.csv",
]


df = pd.read_csv(csv_files[0])
for column in df["CLIMATE"]:
    print("{}: {},".format(column, os.path.isfile(column)))
for column in df["FRF"]:
    print("{}: {},".format(column, os.path.isfile(os.path.join("../front_retreat", column))))

df = pd.read_csv(csv_files[1])
for column in df["CLIMATE"]:
    print("{}: {},".format(column, os.path.isfile(column)))
for column in df["RU"]:
    print("{}: {},".format(column, os.path.isfile(column)))
for column in df["THETA"]:
    print("{}: {},".format(column, os.path.isfile(column)))
