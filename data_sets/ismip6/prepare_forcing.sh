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
ocean_version=v2

grisdir=Projections/GrIS

atmospheredir=Atmosphere_Forcing/aSMB_observed/${climate_version}
atmosphereforcingdir=${basedir}/${grisdir}/${atmospheredir}

oceandir=Ocean_Forcing/Melt_Implementation/${ocean_version}
oceanforcingdir=${basedir}/${grisdir}/${oceandir}

rcm=MARv3.9

start_year=2008
end_year=2100

Climate Forcing
for gcm in MIROC5-rcp26 MIROC5-rcp85 NorESM1-rcp85; do
    for var in aSMB dSMBdz aST dSTdz; do
        eval "cdo -O -f nc4 -z zip_3 mergetime ${atmosphereforcingdir}/${gcm}/${var}/${var}_${rcm}-yearly-${gcm}-{"${start_year}..${end_year}"}.nc ${rcm}_${gcm}-${var}_${start_year}-${end_year}_${climate_version}.nc"
    done
    eval "cdo -O -f nc4 -z zip_3 merge ${rcm}_${gcm}-{aSMB,dSMBdz,aST,dSTdz}_${start_year}-${end_year}_${climate_version}.nc" ${rcm}_${gcm}-climate_${start_year}-${end_year}_${climate_version}.nc""
    adjust_timeline.py -p yearly -a ${start_year}-1-1 -d ${start_year}-1-1  ${rcm}_${gcm}-climate_${start_year}-${end_year}_${climate_version}.nc
    eval "rm ${rcm}_${gcm}-{aSMB,dSMBdz,aST,dSTdz}_${start_year}-${end_year}_${climate_version}.nc"
done

# Ocean Forcing

rcmdir=noresm1-m_rcp8.5
rcmbasename=MAR3.9_NorESM1-M_rcp85

cdo -O -L -f nc4  -z zip_3 merge -chname,basin_runoff,water_input_rate -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${oceanforcingdir}/${rcmdir}/${rcmbasename}_basinRunoff_${ocean_version}.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${oceanforcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing_${ocean_version}.nc ${rcmbasename}_ocean_${ocean_version}.nc
adjust_timeline.py -p yearly -a ${start_year}-1-1 -d ${start_year}-1-1 ${rcmbasename}_ocean_${ocean_version}.nc

rcmdir=miroc-esm-chem_rcp2.6 
rcmbasename=MAR3.9_MIROC-ESM-CHEM_rcp26

cdo -O -L -f nc4  -z zip_3 merge -chname,basin_runoff,water_input_rate -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${oceanforcingdir}/${rcmdir}/${rcmbasename}_basinRunoff_${ocean_version}.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${oceanforcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing_${ocean_version}.nc ${rcmbasename}_ocean_${ocean_version}_${ocean_version}.nc
adjust_timeline.py -p yearly -a ${start_year}-1-1 -d ${start_year}-1-1  ${rcmbasename}_ocean_${ocean_version}.nc

rcmdir=miroc-esm-chem_rcp8.5
rcmbasename=MAR3.9_MIROC-ESM-CHEM_rcp85

cdo -O -L -f nc4  -z zip_3 merge -chname,basin_runoff,water_input_rate -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${oceanforcingdir}/${rcmdir}/${rcmbasename}_basinRunoff_${ocean_version}.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 -setmissval,nan -selyear,${start_year}/${end_year} ${oceanforcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing_${ocean_version}.nc ${rcmbasename}_ocean_${ocean_version}.nc
adjust_timeline.py -p yearly -a ${start_yera}-1-1 -d ${start_yera}-1-1  ${rcmbasename}_ocean_${ocean_version}.nc
