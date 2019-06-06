#!/bin/bash

# Postprocessing initMIP
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

# Split variables
for exp in ctrl asmb; do
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo -f nc4 -z zip_3 splitname  ${odir}_tmp/ex_ismip6_g1000m_v3a_exp_${exp}_2008-1-1_2108-1-1.nc  ${odir}/spatial/
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo  -f nc4 splitname  ${odir}/scalar/ts_ismip6_g1000m_v3a_exp_${exp}_2008-1-1_2108-1-1.nc  ${odir}/${scalar_dir}/
done

# Adjust time axis
for file in ${odir}/spatial/*.nc; do
    ~/pism/util/adjust_timeline.py -i 5 -a ${start_date} -d ${start_date} -p yearly $file
done
for file in ${odir}/${scalar_dir}/*.nc; do
    ~/pism/util/adjust_timeline.py -i 1 -a ${start_date} -d ${start_date} -p yearly $file
done
