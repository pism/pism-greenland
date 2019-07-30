odir=2019_05_vc
grid=900
n=560

for d in gris; do     python historical.py -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_05_hc
grid=900
n=560

for d in gris; do     python historical.py --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

grid=900
n=360

for d in gris; do     python historical.py --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done


# ISMIP 6 Runs

odir=2019_06_calib
n=240

grid=900
for d in gris; do     python historical.py --spatial_ts ismip6 -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

grid=1000
for d in ismip6; do     python historical.py --spatial_ts ismip6 -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

n=360

grid=600
for d in gris; do     python historical.py --spatial_ts ismip6 -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 38:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_06_calib
grid=450
n=560

for d in gris; do  python historical.py -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell  -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_06_calib
n=224

grid=1000
for d in ismip6; do     python historical.py --spatial_ts ismip6 -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_07_hind
n=224

grid=1000
for d in ismip6; do     python historical.py --spatial_ts basic --exstep monthly -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done
