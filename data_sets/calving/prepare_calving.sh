#!/bin/bash


set -x -e

python create_seasonal_calving.py seasonal_calving_1_1980_2020.nc
python create_seasonal_calving.py -s 0 seasonal_calving_0_1980_2020.nc

a=1998-01-01
e=2000-01-01

tct_a=300
tct_e=500

python create_calving_at_thickness.py --threshold_a ${tct_a} --threshold_e ${tct_e} --date_a $a --date_e $e thickness_calving_threshold_${tct_a}_${tct_e}_${a}_${e}.nc
