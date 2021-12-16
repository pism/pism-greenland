#!/bin/bash


set -x -e

python create_seasonal_calving.py seasonal_calving_1_1980_2020.nc
python create_seasonal_calving.py -s 0 seasonal_calving_0_1980_2020.nc
