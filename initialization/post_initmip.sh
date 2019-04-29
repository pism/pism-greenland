#!/bin/bash

# Postprocessing iniMIP
#
# This includes:
# splitname: split files by variable
# settunis: set units to days
# settaxis: set the time axis

odir=$1
IS=GIS
GROUP=UAF
MODEL=PISM

start_date=2008-1-1
scalar_dir=scalar_processed
mkdir -p ${odir}/${scalar_dir}

for exp in ctrl asmb; do
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo -f nc4 -z zip_3 splitname -settunits,day -settaxis,${start_date},00:00:00,5year  ${odir}_tmp/ex_ismip6_g1000m_v3a_exp_${exp}_2008-1-1_2108-1-1.nc  ${odir}/spatial/
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo  -f nc4 splitname -settunits,days  -settaxis,${start_date},0:00:00,1year ${odir}/scalar/ts_ismip6_g1000m_v3a_exp_${exp}_2008-1-1_2108-1-1.nc  ${odir}/${scalar_dir}/
done
