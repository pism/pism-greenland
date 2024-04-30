#!/bin/bash

# Fix grid description, adjust time axis, and set meta data.
#
# The attribute "long\ name" should be "long_name", but I am unable to rename it
# with CDO or NCO.

for file in pr_*.nc; do
    cdo -L -O  setcalendar,365_day -setreftime,1950-1-1 -settaxis,-18050-1-1,0:00:00,50year -setgrid,../../grids/grid_badgeley.txt $file cf_${file}
done

for file in tas_*.nc; do
    cdo -L -O  setattribute,p*@units="Celsius" -setcalendar,365_day -setreftime,1950-1-1 -settaxis,-18050-1-1,0:00:00,50year -setgrid,../../grids/grid_badgeley.txt $file cf_${file}
done
