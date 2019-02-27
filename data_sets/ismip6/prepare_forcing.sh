#!/bin/bash

basedir=$1
grisdir=cryoftp1.gsfc.nasa.gov/ISMIP6/Projections/GrIS
oceandir=Ocean_Forcing/Melt_Implementation
forcingdir=${basedir}/${grisdir}/${oceandir}

rcmdir=miroc-esm-chem_rcp2.6 
rcmbasename=MAR3.9_MIROC-ESM-CHEM_rcp26

cdo -O -L -f nc4  merge -chname,basin_runoff,water_input_rate -setmisstoc,0  ${forcingdir}/${rcmdir}/${rcmbasename}_basinRunoff.nc -chname,thermal_forcing,theta_ocean -setmisstoc,0 ${forcingdir}/${rcmdir}/${rcmbasename}_oceanThermalForcing.nc ${rcmbasename}_ocean.nc
adjust_timeline.py -p yearly -a 1950-1-1 -d 1900-1-1  ${rcmbasename}_ocean.nc
