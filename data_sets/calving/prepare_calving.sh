#!/bin/bash


set -x -e

GRID=3000
infile=../bed_dem/pism_Greenland_${GRID}m_mcb_jpl_v3a_ctrl.nc



# #####################################
# Thickness calving threshold forcing
# #####################################


lat_0=74
lat_1=76

tct_0=0
tct_1=0
outfile=tct_forcing_off.nc
ncks -6 -C -O -v x,y,mask,polar_stereographic $infile $outfile
python tct_forcing.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile
ncatted -a grid_mapping,thickness_calving_threshold,o,c,"polar_stereographic" $outfile

tct_0=200
tct_1=50
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
ncks -6 -C -O -v x,y,mask,polar_stereographic $infile $outfile
python tct_forcing.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile
ncatted -a grid_mapping,thickness_calving_threshold,o,c,"polar_stereographic" $outfile

tct_0=300
tct_1=50
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
ncks -6 -C -O -v x,y,mask,polar_stereographic $infile $outfile
python tct_forcing.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile
ncatted -a grid_mapping,thickness_calving_threshold,o,c,"polar_stereographic" $outfile

tct_0=400
tct_1=50
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
ncks -6 -C -O -v x,y,mask,polar_stereographic $infile $outfile
python tct_forcing.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile
ncatted -a grid_mapping,thickness_calving_threshold,o,c,"polar_stereographic" $outfile

tct_0=500
tct_1=100
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
ncks -6 -C -O -v x,y,mask,polar_stereographic $infile $outfile
python tct_forcing.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile
ncatted -a grid_mapping,thickness_calving_threshold,o,c,"polar_stereographic" $outfile

tct_0=600
tct_1=150
outfile=tct_forcing_${tct_0}myr_${lat_0}n_${tct_1}myr_${lat_1}n.nc
ncks -6 -C -O -v x,y,mask,polar_stereographic $infile $outfile
python tct_forcing.py --tct_0 ${tct_0} --tct_1 ${tct_1} --lat_0 ${lat_0} --lat_1 ${lat_1} $outfile
ncatted -a grid_mapping,thickness_calving_threshold,o,c,"polar_stereographic" $outfile

