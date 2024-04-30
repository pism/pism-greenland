#!/bin/bash

cdo -O -P 6 -f nc4  remapbil,../../grids/gris_ext_g4500m.txt -setgrid,../../grids/bamber_5km_grid.txt -selvar,precipitation,ice_surface_temp pism_Greenland_5km_v1.1.nc tmp.nc
tmp.nc
mpirun -np 6 python ~/pism/bin/fill_missing_petsc.py -v precipitation,ice_surface_temp tmp.nc pism_SeaRISE_SMB_ext_4500m.nc
ncap2 -O -s "air_temp=ice_surface_temp+273.13;"  pism_SeaRISE_SMB_ext_4500m.nc  pism_SeaRISE_SMB_ext_4500m.nc
ncatted -a units,air_temp,o,c,"K" -a units,ice_surface_temp,o,c,"K"  pism_SeaRISE_SMB_ext_4500m.nc
