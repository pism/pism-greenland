#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

import pandas as pd

if __name__ == "__main__":

    T = pd.read_csv("ginr/GINR-S26-Temperature.csv", names=["Year", "Temperature [Celsius]"])
    S = pd.read_csv("ginr/GINR-S26-Salinity.csv", names=["Year", "Salinity [g/kg]"])
    df = pd.concat([T, S]).sort_values(by="Year")
