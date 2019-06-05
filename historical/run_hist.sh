odir=2019_05_vc
grid=900
n=560

for d in gris; do     python historical.py -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_05_hc
grid=900
n=560

for d in gris; do     python historical.py --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_05_hc
grid=4500
n=24

for d in nw; do     python historical.py --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_05_hc
grid=450
n=600

for d in gris; do  python historical.py --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 68:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_05_hc
grid=900
n=240

for d in gris; do     python historical.py --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_05_hcf
grid=450
n=360

for d in jib; do     python historical.py --exstep monthly --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2009-1-1 -q t2standard -s chinook -w 68:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_05_hcpr
grid=450
n=360

for d in jib; do     python historical.py --exstep daily --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2009-1-1 -q t2standard -s chinook -w 68:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done


odir=2019_05_hc_thk
grid=450
n=360

for d in jib; do     python historical.py --exstep daily --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2009-1-1 -q t2standard -s chinook -w 68:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done


grid=900
n=360

for d in gris; do     python historical.py --calving hayhurst_calving -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done


# ISMIP 6 Runs

odir=2019_06_ismip6
n=360

grid=900
for d in gris; do     python historical.py --spatial_ts ismip6 -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

grid=1000
for d in ismip6; do     python historical.py --spatial_ts ismip6 -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

n=360

grid=600
for d in gris; do     python historical.py --spatial_ts ismip6 -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q t2standard -s chinook -w 38:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_routing.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done
