#!/bin/bash

# Create forcing with 4C forcing
ncap2 -O -s "theta_ocean=water_input_rate*0.0 + 4.0;" ../runoff/DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_4500M_DM.nc THETA_1980_2016_EPSG3413_4500M_DM.nc
ncatted -a standard_name,theta_ocean,d,, -a units,theta_ocean,o,c,"Celsius" THETA_1980_2016_EPSG3413_4500M_DM.nc
