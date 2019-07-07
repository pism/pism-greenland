
odir=2019_07_test
n=24



odir=2019_07_test

grid=4500
n=24

for d in gris; do     python variability.py  -d ${d} --o_dir ${odir}  --duration 50000 --step 25000 -q t2small -s chinook -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/variability.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

grid=3600
n=48

for d in gris; do     python variability.py  -d ${d} --o_dir ${odir}  --duration 50000 --step 10000 -q t2small -s chinook -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/variability.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

