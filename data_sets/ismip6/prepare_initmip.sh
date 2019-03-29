#!/bin/bash
# (C) 2019 Andy Aschwanden, University of Alaska Fairbanks
#
# This bash script needs adjust_timeline.py from
# https://github.com/pism/pism
#
# Though the same could probably be done with cdo settimeaxis

set -x 
basedir=$1
climate_version=v2

grisdir=Projections/GrIS

atmospheredir=initMIP/dSMB_epsg3413
atmosphereforcingdir=${basedir}/${grisdir}/${atmospheredir}

asmb_infile=dsmb_01e3413_ISMIP6_${climate_version}.nc

start_year=2008
end_year=2108

# Climate Forcing
rm -rf timline_${start_year}_${end_year}.nc
cdo setattribute,DSMB@units="kg m-2 yr-1" -mulc,910 ${atmosphereforcingdir}/${asmb_infile} ${asmb_infile}
ncrename -O -d x1,x -d y1,y ${asmb_infile} ${asmb_infile}
ncap2 -4 -L 3 -O -s "x=array(-720000,1000,\$x); y=array(-3450000,1000,\$y); defdim(\"time\", 100); time[\$time]=0; for (idx=0;idx < 100; idx++) {time(idx)=idx;};" ${asmb_infile} ${asmb_infile}
ncap2 -4 -L 3 -O -s "ice_surface_temp_anomaly[\$time,\$y,\$x]=0.; climatic_mass_balance_anomaly[\$time,\$y,\$x]=DSMB; for (idx=0;idx < 40; idx++) {climatic_mass_balance_anomaly(idx,:,:)=DSMB(:,:)*(floor(idx) / 40); };"  ${asmb_infile} ${asmb_infile}
adjust_timeline.py -p yearly -a ${start_year}-1-1  -d ${start_year}-1-1  ${asmb_infile}
asmb_outfile=initMIP_climate_forcing_anomalies_${start_year}_${end_year}.nc
ncks -O -4 -L 3 -v DSMB -x  ${asmb_infile}  ${asmb_outfile}
ncatted -a units,ice_surface_temp_anomaly,o,c,"Kelvin"  ${asmb_outfile}
