#!/bin/bash

# Postprocessing prognostic simulations
#
# This includes:
# splitname: split files by variable
# run python postprocessing script

indir=$1
IS=GIS
GROUP=UAF
MODEL=PISM
d=ismip6
g=1000

scalar_dir=scalar_processed
mkdir -p ${indir}/${scalar_dir}

# Split variables
for exp in ctrl_proj 1 5 7 8 9; do
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo -O -f nc4 -z zip_3 splitname  ${indir}_tmp/ex_${d}_g${g}m_v3a_id_EXP-${exp}_2015-1-1_2100-1-1.nc  ${indir}/spatial/
    CDO_FILE_SUFFIX=_${IS}_${GROUP}_${MODEL}_${exp}.nc cdo -O -f nc4 -z zip_3 splitname  ${indir}/scalar/ts_${d}_g${g}m_v3a_id_EXP-${exp}_2015-1-1_2100-1-1.nc  ${indir}/${scalar_dir}/
done

python postprocess_prognostic.py -o ${indir}_processed ${indir}

odir=2019_08_ismip6
mkdir -p $odir/dgmsl
cd $indir/scalar
for file in *.nc; do
    echo $file
    cdo -L setattribute,ice_mass@units="cm" -setattribute,long_mame="contribution to global mean sea level" -divc,362.5 -divc,-1e13 -selvar,limnsw -sub $file -seltimestep,1 $file ../dgmsl/$file
done
