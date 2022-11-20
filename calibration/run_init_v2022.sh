#!/bin/bash

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

    python3 calibrate-v2022.py  --o_dir ${odir} --step 20 --duration 20 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_pseudo_plastic_ctrl.csv ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc 

done

odir=2022_11_tillwat

for (( i=1; i<${nr}+1; i++ )); do
    grid=${grids[$i-1]}
    n=${ns[$i-1]}
    q=${qs[$i-1]}
    w=${ws[$i-1]}

    python3 calibrate-v2022.py  --o_dir ${odir} --step 5 --duration 5 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_pseudo_plastic_ctrl.csv 2022_11_init/state/gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc
done


odir=2022_11_init_ragis

for (( i=1; i<${nr}+1; i++ )); do
    grid=${grids[$i-1]}
    n=${ns[$i-1]}
    q=${qs[$i-1]}
    w=${ws[$i-1]}

    python3 calibrate-v2022.py  --o_dir ${odir} --dataset_version 2022_RAGIS --step 20 --duration 20 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_pseudo_plastic_ctrl.csv ../../best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc 

done


odir=2022_11_tillwat

for (( i=1; i<${nr}+1; i++ )); do
    grid=${grids[$i-1]}
    n=${ns[$i-1]}
    q=${qs[$i-1]}
    w=${ws[$i-1]}

    python3 calibrate-v2022.py  --o_dir ${odir} --step 5 --duration 5 -s chinook -q ${q} -n ${n} -g ${grid} -w ${w} --ensemble_file ../uncertainty_quantification/ensemble_gris_pseudo_plastic_ctrl.csv 2022_11_init_ragis/state/gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc
done

for grid in 1800 1500 1200 900 600 450; do nccopy gris_g${grid}m_v2022_id_CTRL_0_20.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc ; ncks -A ../../../data_sets/velocities/GRE_G${grid}_0000.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc; ncap2 -O -s "where(Band1>100) tillwat=2.0;" gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc; done

for grid in 900; do nccopy gris_g${grid}m_v2022_id_CTRL_0_20.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc ; ncks -A ../../../data_sets/velocities/GRE_G${grid}_0000.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc; ncap2 -O -s "where(Band1>100) tillwat=2.0;" gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc gris_g${grid}m_v2022_id_CTRL_0_20_tillwat.nc; done


for grid in 1800 1500 1200 900 600 450; do nccopy gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc ; ncks -A ../../../data_sets/velocities/GRE_G${grid}_0000.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc; ncap2 -O -s "where(Band1>100) tillwat=2.0;" gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc; done

 for grid in 900; do nccopy gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc ; ncks -A ../../../data_sets/velocities/GRE_G${grid}_0000.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc; ncap2 -O -s "where(Band1>100) tillwat=2.0;" gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc gris_g${grid}m_v2022_RAGIS_id_CTRL_0_20_tillwat.nc; done


