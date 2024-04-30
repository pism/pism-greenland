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

scalar_dir=scalar_processed
mkdir -p ${indir}/${scalar_dir}

# Split variables
for exp in ctrl asmb; do
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo -f nc4 -z zip_3 splitname  ${indir}_tmp/ex_*_exp_${exp}_2008-1-1_2108-1-1.nc  ${indir}/spatial/
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo  -f nc4 -z zip_3 splitname  ${indir}/scalar/ts_*_exp_${exp}_2008-1-1_2108-1-1.nc  ${indir}/${scalar_dir}/
done

python postprocess_init.py -o ${indir}_processed ${indir}
