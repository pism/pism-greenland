#!/bin/bash

ncap2 -O -s "theta_ocean=water_input_rate*0.0 + 277.13;" ../runoff/DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_4500m_DM.nc THETA_1980_2016_EPSG3413_4500m_DM.nc
ncatted -a standard_name,theta_ocean,d,, -a units,theta_ocean,o,c,"K" THETA_1980_2016_EPSG3413_4500m_DM.nc
