#!/bin/bash

set -x

odir=2022_11_tillwat
v=2022
id=CTRL
mkdir -p $odir/flux_gates

for grid in 1800 1500 1200 900; do

python ~/base/pism-analysis/pismanalysis/extract_profiles.py  -v velsurf_mag,uvelsurf,vvelsurf -s  --srs epsg:3413  ~/base/gris-analysis/flux-gates/greenland-flux-gates-29-250m.shp ${odir}/state/gris_g${grid}m_v${v}_id_${id}_0_5.nc ${odir}/flux_gates/flux_gate_29_250m_gris_g${grid}m_v${v}_id_${id}_0_5.nc

ncap2 -O -s "velsurf_normal=(uvelsurf*nx + vvelsurf*ny); config[]=\" \"; config@grid_dx_meters=${grid}; config@init=\"TILLWAT\"; config@bed=\"v5\";" ${odir}/flux_gates/flux_gate_29_250m_gris_g${grid}m_v${v}_id_${id}_0_5.nc ${odir}/flux_gates/flux_gate_29_250m_gris_g${grid}m_v${v}_id_${id}_0_5.nc

done


odir=2022_11_init
v=2022
id=CTRL
mkdir -p $odir/flux_gates

for grid in 1800 1500 1200 900; do

python ~/base/pism-analysis/pismanalysis/extract_profiles.py  -v velsurf_mag,uvelsurf,vvelsurf -s  --srs epsg:3413  ~/base/gris-analysis/flux-gates/greenland-flux-gates-29-250m.shp ${odir}/state/gris_g${grid}m_v${v}_id_${id}_0_20.nc ${odir}/flux_gates/flux_gate_29_250m_gris_g${grid}m_v${v}_id_${id}_0_20.nc

ncap2 -O -s "velsurf_normal=(uvelsurf*nx + vvelsurf*ny); config[]=\" \"; config@grid_dx_meters=${grid}; config@init=\"CTRL\"; config@bed=\"v5\";" ${odir}/flux_gates/flux_gate_29_250m_gris_g${grid}m_v${v}_id_${id}_0_20.nc ${odir}/flux_gates/flux_gate_29_250m_gris_g${grid}m_v${v}_id_${id}_0_20.nc

done


odir=best_v1

for grid in 4500 3600 1800 1500 1200 900; do

python ~/base/pism-analysis/pismanalysis/extract_profiles.py -v velsurf_mag,uvelsurf,vvelsurf -s --srs epsg:3413  ~/base/gris-analysis/flux-gates/greenland-flux-gates-29-250m.shp /Volumes/zachariae/best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc ../data_sets/flux_gates/flux_gate_29_250m_g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc

ncap2 -O -s "velsurf_normal=(uvelsurf*nx + vvelsurf*ny); config[]=\" \"; config@grid_dx_meters=${grid}; config@init=\"CTRL\"; config@bed=\"v1\";" ../data_sets/flux_gates/flux_gate_29_250m_g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc ../data_sets/flux_gates/flux_gate_29_250m_g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc
done


for grid in 600; do

python ~/base/pism-analysis/pismanalysis/extract_profiles.py -v velsurf_mag,uvelsurf,vvelsurf -s --srs epsg:3413  ~/base/gris-analysis/flux-gates/greenland-flux-gates-29-250m.shp /Volumes/zachariae/best_v1/g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc ../data_sets/flux_gates/flux_gate_29_250m_g${grid}m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc

done

python ~/base/pism-analysis/pismanalysis/flux_gate_analysis.py -v velsurf_normal --label_params bed --obs_file ../data_sets/velocities/flux_gates_29_250m_GRE_G0120_0000.nc ../data_sets/flux_gates/flux_gate_29_250m_g900m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc 2022_11_tillwat/flux_gates/flux_gate_29_250m_gris_g900m_v2022_id_CTRL_0_5.nc
