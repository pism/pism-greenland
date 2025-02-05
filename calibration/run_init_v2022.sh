#!/bin/bash

grid=1800
 for v in {1..5}; do PISM_PREFIX=~/local/pism-dev/ python calibrate-v2022.py --spatial_ts none --o_dir 2024_05_init_bedmachine --dataset_version ${v}_GrIMP --b wc --duration 20 -s pleiades_haswell -q long -n 120 -g $grid -w 24:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv  ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc; qsub /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine/run_scripts/cc_g${grid}m_v${v}_GrIMP_id_BAYES-MEDIAN_0_20.sh; qsub /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine/run_scripts/cc_g${grid}m_v${v}_GrIMP_id_CTRL_0_20.sh; done

grid=900
 for v in {1..5}; do PISM_PREFIX=~/local/pism-dev/ python calibrate-v2022.py --spatial_ts none --o_dir 2024_05_init_bedmachine --dataset_version ${v}_GrIMP --b wc --duration 20 -s pleiades_haswell -q long -n 240 -g $grid -w 24:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv  ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc; qsub /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine/run_scripts/cc_g${grid}m_v${v}_GrIMP_id_BAYES-MEDIAN_0_20.sh; qsub /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine/run_scripts/cc_g${grid}m_v${v}_GrIMP_id_CTRL_0_20.sh; done


grid=900
for e in CTRL BAYES-MEDIAN; do
for v in {1..5}; do PISM_PREFIX=~/local/pism-dev/ python calibrate-v2022.py --spatial_ts none --o_dir 2024_05_init_bedmachine_tw --dataset_version ${v}_GrIMP --b wc --duration 20 -s pleiades_haswell -q long -n 240 -g $grid -w 16:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv   /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine/state/gris_g${grid}m_v${v}_GrIMP_id_${e}_0_20_tillwat.nc; qsub /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine_tw/run_scripts/cc_g${grid}m_v${v}_GrIMP_id_${e}_0_20.sh;
done
done

grid=1800
for e in CTRL BAYES-MEDIAN; do
for v in {1..5}; do PISM_PREFIX=~/local/pism-dev/ python calibrate-v2022.py --spatial_ts none --o_dir 2024_05_init_bedmachine_tw --dataset_version ${v}_GrIMP --b wc --duration 20 -s pleiades_haswell -q long -n 120 -g $grid -w 12:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv   /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine/state/gris_g${grid}m_v${v}_GrIMP_id_${e}_0_20_tillwat.nc; qsub /nobackupp17/aaschwan/pism-greenland/calibration/2024_05_init_bedmachine_tw/run_scripts/cc_g${grid}m_v${v}_GrIMP_id_${e}_0_20.sh;
done
done




PISM_PREFIX=~/local/pism-dev python calibrate-v2022.py --spatial_ts basic --o_dir 2024_04_init_krig_tw --dataset_version 2024_KRIG --b wc --duration 20 -s pleiades_haswell -n 240 -q long  -g 900 -w 24:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv  2024_04_init_krig/state/gris_g900m_v2024_KRIG_id_CTRL_0_20_tillwat.nc



PISM_PREFIX=$HOME/local-rl8/pism python calibrate-v2022.py --o_dir 2023_10_init_grimp --dataset_version 2023_GIMP --b wc --duration 25 -s chinook-rl8 -q t2standard -n 200 -g 600 -w 48:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv ../../best_v1/g600m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc

PISM_PREFIX=$HOME/local-rl8/pism python calibrate-v2022.py --o_dir 2023_10_init_ragis --dataset_version 2023_RAGIS --b wc --duration 25 -s chinook-rl8 -q t2standard -n 200 -g 600 -w 48:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv ../../best_v1/g600m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc
i


grid=1200
PISM_PREFIX=$HOME/local-rl8/pism/ python calibrate-v2022.py --o_dir 2023_09_pseudo_plastic --dataset_version 2023_GIMP --b wc --spatial_ts basic --exstep 1 --step 25 --duration 25 -s chinook-rl8 -q t2small -n 80 -g $grid -w 48:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_speed_calibration_pseudo_plastic_lhs_1000.csv ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc


grid=1200
PISM_PREFIX=$HOME/local-rl8/pism/ python calibrate-v2022.py --o_dir 2023_09_pseudo_plastic_test --dataset_version 2023_GIMP --b wc --spatial_ts basic --exstep 1 --step 25 --duration 25 -s chinook-rl8 -q t2small -n 80 -g $grid -w 36:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_speed_calibration_pseudo_plastic_lhs_1000.csv  ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc

grid=1200
PISM_PREFIX=$HOME/local-rl8/pism/ python calibrate-v2022.py --o_dir 2023_09_pseudo_plastic_test --dataset_version 2023_GIMP --b wc --spatial_ts basic --exstep 1 --step 25 --duration 25 -s chinook-rl8 -q debug -n 40 -g $grid -w 36:00:00 --ensemble_file ../uncertainty_quantification/ensemble_gris_speed_calibration_pseudo_plastic_lhs_1000.csv  ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc

# Two steps:
# 1. Update to v2022 bed/ice thickness
# 2. Update tillwat

declare -a grids=(1800 1500 1200 900 600 450)
declare -a ns=(24 24 48 72 96 144)
declare -a qs=("t2small" "t2small" "t2small" "t2standard" "t2standard" "t2standard")
declare -a ws=("2:00:00" "4:00:00" "4:00:00" "6:00:00" "28:00:00" "72:00:00")
nr=${#grids[@]}

odir=2022_11_init

for (( i=1; i<${nr}+1; i++ )); do
    grid=${grids[$i-1]}
    n=${ns[$i-1]}
    q=${qs[$i-1]}
    w=${ws[$i-1]}

    python3 calibrate-v2022.py  --o_dir ${odir} --step 20 --duration 20 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc 

done

odir=2022_11_tillwat

for (( i=1; i<${nr}+1; i++ )); do
    grid=${grids[$i-1]}
    n=${ns[$i-1]}
    q=${qs[$i-1]}
    w=${ws[$i-1]}

    python3 calibrate-v2022.py  --o_dir ${odir} --step 5 --duration 5 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv 2022_11_init/state/gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc
done


odir=2023_02_init_ragis

for (( i=1; i<${nr}+1; i++ )); do
    grid=${grids[$i-1]}
    n=${ns[$i-1]}
    q=${qs[$i-1]}
    w=${ws[$i-1]}

    python3 calibrate-v2022.py  --o_dir ${odir} --dataset_version 2023_RAGIS --step 20 --duration 20 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc 

done


odir=2022_11_tillwat

for (( i=1; i<${nr}+1; i++ )); do
    grid=${grids[$i-1]}
    n=${ns[$i-1]}
    q=${qs[$i-1]}
    w=${ws[$i-1]}

    python3 calibrate-v2022.py  --o_dir ${odir} --step 5 --duration 5 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_ctrl.csv 2022_11_init_ragis/state/gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc
done

for grid in 1800 1500 1200 900 600 450; do nccopy gris_g${grid}m_v2022_id_CTRL_0_20.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc ; ncks -A ../../../data_sets/velocities/GRE_G${grid}_0000.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc; ncap2 -O -s "where(Band1>100) tillwat=2.0;" gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc; done


for grid in 1800 1500 1200 900 600 450; do nccopy gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc ; ncks -A ../../../data_sets/velocities/GRE_G${grid}_0000.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc; ncap2 -O -s "where(Band1>100) tillwat=2.0;" gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc; done


grid=1800
n=68
s=stampede2
q=normal
w=4:00:00
dataset=2022
uq=pseudo_plastic_lhs_500
odir=$SCRATCH/2022_11_${uq}_test

PISM_PREFIX=$HOME/local/pism-dev/bin python3 calibrate-v2022.py  --o_dir ${odir} --dataset_version $dataset --step 20 --duration 20 -s $s -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_${uq}.csv ../data_sets/initial_states/gris_g${grid}m_v${dataset}_id_CTRL_0_5.nc

grid=1800
n=24
s=stampede2
q=normal
w=4:00:00
dataset=2022
uq=pseudo_plastic_lhs_500
odir=$SCRATCH/2022_11_${uq}_test_n24

PISM_PREFIX=$HOME/local/pism-dev/bin python3 calibrate-v2022.py  --o_dir ${odir} --dataset_version $dataset --step 1 --duration 1 -s $s -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_${uq}.csv ../data_sets/initial_states/gris_g${grid}m_v${dataset}_id_CTRL_0_5.nc

grid=1800
n=24
s=stampede2
q=normal
w=4:00:00
dataset=2022
uq=pseudo_plastic_lhs_500
odir=$SCRATCH/2022_11_${uq}_test_n24_icc

PISM_PREFIX=$HOME/local/pism-dev/bin python3 calibrate-v2022.py  --o_dir ${odir} --dataset_version $dataset --step 1 --duration 1 -s $s -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_${uq}.csv ../data_sets/initial_states/gris_g${grid}m_v${dataset}_id_CTRL_0_5.nc

grid=1500
n=68
s=stampede2
q=normal
w=4:00:00
dataset=2022
uq=pseudo_plastic_lhs_500
odir=$SCRATCH/2022_11_${uq}_test

PISM_PREFIX=$HOME/local/pism-dev/bin python3 calibrate-v2022.py  --o_dir ${odir} --dataset_version $dataset --step 20 --duration 20 -s $s -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_${uq}.csv ../data_sets/initial_states/gris_g${grid}m_v${dataset}_id_CTRL_0_5.nc

grid=1800
n=24
s=chinook
q=t2small
w=4:00:00
dataset=2022
uq=pseudo_plastic_lhs_500
odir=2022_11_${uq}_test

PISM_PREFIX=$HOME/local/pism/bin python3 calibrate-v2022.py  --o_dir ${odir} --dataset_version $dataset --step 1 --duration 1 -s $s -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_${uq}.csv ../data_sets/initial_states/gris_g${grid}m_v${dataset}_id_CTRL_0_5.nc
