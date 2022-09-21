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
        REMAP_EXTRAPOLATE=on cdo -L -P $N -O setmisstoc,1  -remapbil,../../grids/gris_ext_g4500m.txt  -setattribute,frac_P@units="1" -chname,posterior_${pctl},frac_P -selvar,posterior_${pctl} -settbounds,50year -setcalendar,365_day -setreftime,1-1-1 -settaxis,-19999-1-1,0:00:00,50year -setgrid,../../grids/grid_badgeley.txt pr_${r}_Badgeley_etal_2020.nc pr_Badgeley_etal_2020_id_${r}-${pctl}.nc
        # CDO sets "1" units as a numeric value, and PISM bails...
        ncatted -a units,time,o,c,"365 days since 1-1-1 00:00:00" -a units,frac_P,o,c,"1" pr_Badgeley_etal_2020_id_${r}-${pctl}.nc
        ncap2 -O -s "*sz_idx=time.size(); for(*idx=0;idx<sz_idx;idx++){time_bnds(idx,0)=time(idx)-50; time_bnds(idx,1)=time(idx);};" pr_Badgeley_etal_2020_id_${r}-${pctl}.nc pr_Badgeley_etal_2020_id_${r}-${pctl}.nc
    done
    for r in main S1 S2 S3 S4; do
        REMAP_EXTRAPOLATE=on cdo -L -P $N -O setmisstoc,0 -remapbil,../../grids/gris_ext_g4500m.txt -setattribute,delta_T@units="K" -chname,posterior_${pctl},delta_T -selvar,posterior_${pctl} -settbounds,50year -setcalendar,365_day -setreftime,1-1-1 -settaxis,-19999-1-1,0:00:00,50year -setgrid,../../grids/grid_badgeley.txt tas_${r}_Badgeley_etal_2020.nc tas_Badgeley_etal_2020_id_${r}-${pctl}.nc
        ncap2 -O -s "*sz_idx=time.size(); for(*idx=0;idx<sz_idx;idx++){time_bnds(idx,0)=time(idx)-50; time_bnds(idx,1)=time(idx);};" tas_Badgeley_etal_2020_id_${r}-${pctl}.nc tas_Badgeley_etal_2020_id_${r}-${pctl}.nc
        cdo -O --reduce_dim fldmean tas_Badgeley_etal_2020_id_${r}-${pctl}.nc tas_Badgeley_etal_2020_id_${r}-${pctl}_fldmean.nc
        ncatted -a units,time,o,c,"365 days since 1-1-1 00:00:00" tas_Badgeley_etal_2020_id_${r}-${pctl}.nc
        ncatted -a units,time,o,c,"365 days since 1-1-1 00:00:00" tas_Badgeley_etal_2020_id_${r}-${pctl}_fldmean.nc
    done
done

            
