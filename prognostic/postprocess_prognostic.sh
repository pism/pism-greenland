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

declare -a oldexps=("ctrl_proj" "1" "2" "3" "4" "5" "6" "7" "8" "9" "A1" "A2" "A3" "A4" "B1" "B2" "B3" "B4" "B5" "C1" "C2" "C3" "C4" "C5" "C6" "C9" "C8" "C9" "C10" "D1" "D2" "D3" "D4" "D5" "D6" "D7" "D8" "D9" "D10" "D11" "D12" "D13" "D14" "D15" "D16" "D17" "D18" "D19" "D20" "D21" "D22")
declare -a newexps=("ctrl_proj" "exp01" "exp02" "exp03" "exp04" "exp05" "exp06" "exp07" "exp08" "exp09" "expa01" "expa02" "expa03" "expa04" "expb01" "expb02" "expb03" "expb04" "expb05" "expc01" "expc02" "expc03" "expc04" "expc05" "expc06" "expc07" "expc08" "expc09" "expc10" "expd01" "expd02" "expd03" "expd04" "expd05" "expd06" "expd07" "expd08" "expd09" "expd10"  "expd11" "expd12" "expd13" "expd14" "expd15" "expd16" "expd17" "expd18" "expd19" "expd20" "expd21" "expd22")
n=${#oldexps[@]}



for (( i=1; i<${n}+1; i++ )); do
    oldexp=${oldexps[$i-1]}
    newexp=${newexps[$i-1]}
    echo "Renaming ${oldexp} to ${newexp}"
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${newexp}.nc cdo -O -f nc4 -z zip_3 splitname  ${indir}_tmp/ex_${d}_g${g}m_v3a_id_EXP-${oldexp}_2015-1-1_2100-1-1.nc  ${indir}/spatial/
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${newexp}.nc cdo -O -f nc4 -z zip_3 splitname  ${indir}/scalar/ts_${d}_g${g}m_v3a_id_EXP-${oldexp}_2015-1-1_2100-1-1.nc  ${indir}/${scalar_dir}/
done

# python postprocess_prognostic.py --model $MODEL -o ${outdir} ${indir}
