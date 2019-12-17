#!/bin/bash

# Postprocessing prognostic simulations
#
# This includes:
# splitname: split files by variable
# run python postprocessing script

set -x

indir=$1
outdir=$2
IS=GIS
GROUP=UAF
MODEL=$3
d=ismip6
g=1000

scalar_dir=scalar_processed
mkdir -p ${indir}/${scalar_dir}

declare -a oldexps=("CALIB-G1000M")
declare -a newexps=("historical")
n=${#oldexps[@]}



for (( i=1; i<${n}+1; i++ )); do
    oldexp=${oldexps[$i-1]}
    newexp=${newexps[$i-1]}
    echo "Renaming ${oldexp} to ${newexp}"
    mkdir -p ${outdir}/${GROUP}/${MODEL}/${newexp}_01
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${newexp}.nc cdo -O -f nc4 -z zip_3 splitname  ${indir}_tmp/ex_${d}_g${g}m_v3a_id_${oldexp}_2008-1-1_2015-1-1.nc  ${indir}/spatial/
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${newexp}.nc cdo -O -f nc4 -z zip_3 splitname  ${indir}/scalar/ts_${d}_g${g}m_v3a_id_${oldexp}_2008-1-1_2015-1-1.nc  ${indir}/${scalar_dir}/
done

python ../prognostic/postprocess_prognostic.py --model $MODEL -o ${outdir} ${indir}
