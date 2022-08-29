#!/bin/bash

# Fix grid description, adjust time axis, and set meta data.
#
# The attribute "long\ name" should be "long_name", but I am unable to rename it
# with CDO or NCO.
set -x
set -e

N=4

for pctl in 5 mean 95; do
    for r in low main high; do
        REMAP_EXTRAPOLATE=on cdo -P $N -O setmisstoc,1  -remapbil,../../grids/gris_ext_g4500m.txt  -setattribute,units@frac_P="1" -chname,posterior_${pctl},frac_P -selvar,posterior_${pctl} -setcalendar,365_day -setreftime,1-1-1 -settaxis,-19999-1-1,0:00:00,50year -setgrid,../../grids/grid_badgeley.txt pr_${r}_Badgeley_etal_2020.nc pr_Badgeley_etal_2020_id_${r}-${pctl}.nc
    done
    for r in main S1 S2 S3 S4; do
        REMAP_EXTRAPOLATE=on cdo -P $N -O setmisstoc,0 -remapbil,../../grids/gris_ext_g4500m.txt -setattribute,units@delta_T="K" -chname,posterior_${pctl},delta_T -selvar,posterior_${pctl} -setcalendar,365_day -setreftime,1-1-1 -settaxis,-19999-1-1,0:00:00,50year -setgrid,../../grids/grid_badgeley.txt tas_${r}_Badgeley_etal_2020.nc tas_Badgeley_etal_2020_id_${r}-${pctl}.nc
        cdo -O fldmean tas_Badgeley_etal_2020_id_${r}-${pctl}.nc tas_Badgeley_etal_2020_id_${r}-${pctl}_fldmean.nc
    done
done

            
