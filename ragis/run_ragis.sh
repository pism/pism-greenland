n=24
grid=1800
odir=2022_08_init

python hindcast.py -e ../uncertainty_quantification/ensemble_gris_ragis_calving_lhs_20_w_posterior.csv  -w 12:00:00 -s chinook -q t2small -g ${grid} -n${n} --start 1980-1-1 --end 1990-1-1 --o_dir ${odir} ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc
