#!/bin/bash
# (C) 2019 Andy Aschwanden, University of Alaska Fairbanks
#
# This bash script needs adjust_timeline.py from
# https://github.com/pism/pism
#
# Though the same could probably be done with cdo settimeaxis

set -x 
basedir=$1
climate_version=v1
ocean_version=v3

grisdir=Projections/GrIS

atmospheredir=Atmosphere_Forcing/aSMB_observed/${climate_version}
atmosphereforcingdir=${basedir}/${grisdir}/${atmospheredir}

oceandir=Ocean_Forcing/Melt_Implementation/${ocean_version}
oceanforcingdir=${basedir}/${grisdir}/${oceandir}

rcm=MARv3.9

start_year=1990
end_year=2100

# # Climate Forcing
for gcm in MIROC5-rcp26 MIROC5-rcp85 NorESM1-rcp85; do
    for var in aSMB dSMBdz aST dSTdz; do
        eval "cdo -O -f nc4 -z zip_3 mergetime ${atmosphereforcingdir}/${gcm}/${var}/${var}_${rcm}-yearly-${gcm}-{"${start_year}..${end_year}"}.nc ${rcm}_${gcm}-${var}_${start_year}-${end_year}_${climate_version}.nc"
    done
    eval "cdo -O -f nc4 -z zip_3 chname,aSMB,climatic_mass_balance_anomaly,dSMBdz,climatic_mass_balance_gradient,aST,ice_surface_temp_anomaly,dSTdz,ice_surface_temp_gradient -merge ${rcm}_${gcm}-{aSMB,dSMBdz,aST,dSTdz}_${start_year}-${end_year}_${climate_version}.nc" ${rcm}_${gcm}_climate_${start_year}-${end_year}_${climate_version}.nc""
    adjust_timeline.py -p yearly -a ${start_year}-1-1 -d ${start_year}-1-1  ${rcm}_${gcm}_climate_${start_year}-${end_year}_${climate_version}.nc
    eval "rm ${rcm}_${gcm}-{aSMB,dSMBdz,aST,dSTdz}_${start_year}-${end_year}_${climate_version}.nc"
done

# Ocean Forcing

declare -a rcmdirs=("csiro-mk3.6_rcp8.5" "hadgem2-es_rcp8.5" "ipsl-cm5-mr_rcp8.5" "miroc-esm-chem_rcp2.6" "miroc-esm-chem_rcp8.5" "noresm1-m_rcp8.5")
declare -a rcmbasenames=("MAR3.9_CSIRO-Mk3.6_rcp85" "MAR3.9_HadGEM2-ES_rcp85" "MAR3.9_IPSL-CM5-MR_rcp85" "MAR3.9_MIROC-ESM-CHEM_rcp26" "MAR3.9_MIROC-ESM-CHEM_rcp85" "MAR3.9_NorESM1-M_rcp85")
n=${#rcmdirs[@]}

for (( i=1; i<${n}+1; i++ )); do
    rcmdir=${rcmdirs[$i-1]}
    rcmbasename=${rcmbasenames[$i-1]}
    nccopy ${oceanforcingdir}/${rcmdir}/${rcmbasename}_basinRunoff_${ocean_version}.nc ${rcmbasename}_basinRunoff_${ocean_version}.nc
    nccopy ${oceanforcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing_${ocean_version}.nc ${rcmbasename}_oceanThermalForcing_${ocean_version}.nc
    cdo -O -f nc4  -z zip_3 merge -chname,basin_runoff,water_input_rate -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${rcmbasename}_basinRunoff_${ocean_version}.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${rcmbasename}_oceanThermalForcing_${ocean_version}.nc ${rcmbasename}_ocean_${start_year}-${end_year}_${ocean_version}.nc
    adjust_timeline.py -p yearly -a ${start_year}-1-1 -d ${start_year}-1-1  ${rcmbasename}_ocean_${start_year}-${end_year}_${ocean_version}.nc
done
