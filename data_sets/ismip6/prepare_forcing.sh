#!/bin/bash

set -x 
basedir=$1
version=v1
grisdir=Projections/GrIS
atmospheredir=Atmosphere_Forcing/aSMB_observed/${version}/
oceandir=Ocean_Forcing/Melt_Implementation
atmosphereforcingdir=${basedir}/${grisdir}/${atmospheredir}
oceanforcingdir=${basedir}/${grisdir}/${oceandir}
rcm=MARv3.9

start_year=1980
end_year=1982

# for gcm in MIROC5-rcp26 MIROC5-rcp85 NorESM1-rcp85; do
for gcm in MIROC5-rcp26; do
    for var in aSMB dSMBdz aST dSTdz; do
        eval "cdo -O -f nc4 -z zip_3 mergetime ${atmosphereforcingdir}/${gcm}/${var}/${var}_${rcm}-yearly-${gcm}-{"${start_year}..${end_year}"}.nc ${rcm}_${gcm}-${var}_${start_year}-${end_year}.nc"
    done
    eval "cdo -O -f nc4 -z zip_3 merge ${rcm}_${gcm}-{aSMB,dSMBdz,aST,dSTdz}_${start_year}-${end_year}.nc" ${rcm}_${gcm}-climate_${start_year}-${end_year}.nc""
    eval "rm ${rcm}_${gcm}-{aSMB,dSMBdz,aST,dSTdz}_${start_year}-${end_year}.nc"
done





rcmdir=noresm1-m_rcp8.5
rcmbasename=MAR3.9_NorESM1-M_rcp85

cdo -O -L -f nc4  -z zip_3 merge -chname,basin_runoff,water_input_rate -setmisstoc,0 -setmissval,nan ${oceanforcingdir}/${rcmdir}/${rcmbasename}_basinRunoff.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 -setmissval,nan ${oceanforcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing.nc ${rcmbasename}_ocean.nc
adjust_timeline.py -p yearly -a 1950-1-1 -d 1900-1-1  ${rcmbasename}_ocean.nc

rcmdir=miroc-esm-chem_rcp2.6 
rcmbasename=MAR3.9_MIROC-ESM-CHEM_rcp26

cdo -O -L -f nc4  -z zip_3 merge -chname,basin_runoff,water_input_rate -setmisstoc,0 -setmissval,nan ${oceanforcingdir}/${rcmdir}/${rcmbasename}_basinRunoff.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 -setmissval,nan ${oceanforcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing.nc ${rcmbasename}_ocean.nc
adjust_timeline.py -p yearly -a 1950-1-1 -d 1900-1-1  ${rcmbasename}_ocean.nc

rcmdir=miroc-esm-chem_rcp8.5
rcmbasename=MAR3.9_MIROC-ESM-CHEM_rcp85

cdo -O -L -f nc4  -z zip_3 merge -chname,basin_runoff,water_input_rate -setmisstoc,0 -setmissval,nan ${oceanforcingdir}/${rcmdir}/${rcmbasename}_basinRunoff.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 -setmissval,nan ${oceanforcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing.nc ${rcmbasename}_ocean.nc
adjust_timeline.py -p yearly -a 1950-1-1 -d 1900-1-1  ${rcmbasename}_ocean.nc
