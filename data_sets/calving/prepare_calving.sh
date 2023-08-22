#!/bin/bash


set -x -e

for h in 1.25 1.5 3.0; do
    python create_seasonal_calving.py --calving_high ${h} seasonal_calving_${h}_2001_1980_2020.nc
done

# #####################################
# Thickness calving threshold forcing
# #####################################


lat_0=74
lat_1=76

tct_0=200
tct_1=50
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
python create_calving_at_thickness.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile

tct_0=300
tct_1=50
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
python create_calving_at_thickness.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile

tct_0=400
tct_1=50
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
python create_calving_at_thickness.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile


