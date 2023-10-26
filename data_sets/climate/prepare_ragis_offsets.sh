#!/bin/bash

rm  ragis_dem_climate_1980_2020.nc
create_timeline.py -a 1980-01-01 -e 2020-01-01 -d 1980-01-01 ragis_dem_climate_1980_2020.nc
ncap2 -O -s "delta_T[\$time]=-2.0f; frac_P[\$time]=1.2" ragis_dem_climate_1980_2020.nc ragis_dem_climate_1980_2020.nc
ncatted -a units,delta_T,o,c,"K" -a units,frac_P,o,c,"" ragis_dem_climate_1980_2020.nc
