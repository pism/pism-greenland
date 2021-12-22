#!/bin/bash

odir=2021_12_ctrl
grid=600
ensfile=jib_ctrl.csv

for file in jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.00_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.25_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.50_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_SEASON-CALV-0.4-200-1.50_1980-1-1_2010-1-1; do
    qsub postprocess_jib.sh $odir $file $grid $ensfile
done

mkdir -p $odir/csv

python  ../uncertainty_quantification/nc2csv.py -o $odir/csv/fldmean_ts.csv $odir/processed/fldmean_masked_ex_jib_g600m_v1_RAGIS_id_*
python  ../uncertainty_quantification/nc2csv.py -o $odir/csv/fldsum_ts.csv $odir/processed/fldsum_masked_ex_jib_g600m_v1_RAGIS_id_*


odir=2021_12_ctrl
grid=600
ensfile=jib_ctrl.csv

for file in jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.00_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.25_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.50_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_SEASON-CALV-0.4-200-1.50_1980-1-1_2010-1-1; do
    qsub postprocess_jib.sh $odir $file $grid $ensfile
done

mkdir -p $odir/csv

python  ../uncertainty_quantification/nc2csv.py -o $odir/csv/fldmean_ts.csv $odir/processed/fldmean_masked_ex_jib_g600m_v1_RAGIS_id_*
python  ../uncertainty_quantification/nc2csv.py -o $odir/csv/fldsum_ts.csv $odir/processed/fldsum_masked_ex_jib_g600m_v1_RAGIS_id_*


odir=2021_12_ctrl
grid=600
ensfile=jib_ctrl.csv

for file in jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.00_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.25_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_CTRL-0.4-200-1.50_1980-1-1_2010-1-1 jib_g600m_v1_RAGIS_id_SEASON-CALV-0.4-200-1.50_1980-1-1_2010-1-1; do
    qsub postprocess_jib.sh $odir $file $grid $ensfile
done

mkdir -p $odir/csv

python  ../uncertainty_quantification/nc2csv.py -o $odir/csv/fldmean_ts.csv $odir/processed/fldmean_masked_ex_jib_g600m_v1_RAGIS_id_*
python  ../uncertainty_quantification/nc2csv.py -o $odir/csv/fldsum_ts.csv $odir/processed/fldsum_masked_ex_jib_g600m_v1_RAGIS_id_*

odir=2021_12_fractures
grid=600
ensfile=jib_fractures.csv

for id in {0..79}; do
    file=jib_g600m_v1_RAGIS_id_${id}_1980-1-1_2010-1-1
    if [ -f "${odir}/state/$file" ]; then
        echo $file
        qsub postprocess_jib.sh $odir $file $grid $ensfile
    fi
done

mkdir $odir/processed_scalar
for id in {0..79}; do
    file=ts_jib_g600m_v1_RAGIS_id_${id}_1980-1-1_2010-1-1
    if [ -f "${odir}/scalar/${file}.nc" ]; then
        cdo -O yearmean $odir/scalar/${file}.nc $odir/processed_scalar/${file}_YM.nc
    fi
done

mkdir -p $odir/csv

python  ../util/nc2csv.py -o $odir/csv/fldmean_ts.csv $odir/processed/fldmean_masked_ex_jib_g600m_v1_RAGIS_id_*
python  ../util/nc2csv.py -o $odir/csv/fldsum_ts.csv $odir/processed/fldsum_masked_ex_jib_g600m_v1_RAGIS_id_*


python  ../util/nc2csv.py -o $odir/csv/scalar_ts.csv $odir/processed_scalar/ts_jib_g600m_v1_RAGIS_id_*.nc
