#!/bin/bash

# Postprocessing initMIP
#
# This includes:
# splitname: split files by variable
# settunis: set units to days
# settaxis: set the time axis

indir=$1
IS=GIS
GROUP=UAF
MODEL=PISM
d=ismip6
g=1000

scalar_dir=scalar_processed
mkdir -p ${indir}/${scalar_dir}

# Split variables
for exp in 1; do
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo -f nc4 -z zip_3 splitname  ${indir}_tmp/ex_${d}_g${g}m_v3a_id_EXP-${exp}-VCM-CALIB-G1000M_2015-1-1_2100-1-1.nc  ${indir}/spatial/
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo  -f nc4 -z zip_3 splitname  ${indir}/scalar/ts_${d}_g${g}m_v3a_id_EXP-${exp}-VCM-CALIB-G1000M_2015-1-1_2100-1-1.nc  ${indir}/${scalar_dir}/
done

# python postprocess_prognostic.py -o ${indir}_processed ${indir}
