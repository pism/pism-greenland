
uq=fractures_all
odir=2022_06_${uq}
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2020-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_02_init_2/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.6-300-1.5_1980-1-1_1990-1-1.nc


uq=full
odir=2022_06_${uq}
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1984-1-1 -q t2small -s chinook -w 6:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_02_init_2/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.6-300-1.5_1980-1-1_1990-1-1.nc

uq=fractures_all
odir=2022_06_${uq}_float
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 2000-1-1 --end 2020-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_06_fractures_all/state/jib_g600m_v1_RAGIS_id_14_1980-1-1_2000-1-1.nc


odir=2021_12_init
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/jib_init.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc

odir=2022_01_init
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/jib_init.csv 2021_12_init/state/jib_g${grid}m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc

odir=2022_02_init
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/jib_init.csv 2021_12_init/state/jib_g${grid}m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc

odir=2022_02_init_2
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/jib_init.csv 2022_02_init/state/jib_g${grid}m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc

uq=fractures
odir=2022_03_${uq}
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2020-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_02_init_2/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.6-100-1.5_1980-1-1_1990-1-1.nc 


uq=fractures_all
odir=2022_06_${uq}
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2020-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_02_init_2/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.6-300-1.5_1980-1-1_1990-1-1.nc


odir=2021_12_init_blatter
grid=600
d=jib
n=24

python historical.py --stress_balance blatter --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 68:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/jib_init.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc

odir=2021_12_init
grid=450
d=jib
n=48

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 18:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/jib_init.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc


uq=fractures_steady
odir=2022_01_${uq}
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_01_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.5-250-1.5_1980-1-1_1990-1-1.nc

uq=fractures_steady
odir=2022_01_${uq}
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_01_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.5-250-1.5_1980-1-1_1990-1-1.nc

uq=calving
odir=2022_01_${uq}
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_${uq}.csv 2022_01_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.5-250-1.5_1980-1-1_1990-1-1.nc


odir=2021_12_ctrl
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/jib_ctrl.csv 2021_12_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc

odir=2021_12_fractures
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_fractures.csv 2021_12_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc

odir=2021_12_fractures
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_fractures.csv 2021_12_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc

odir=2021_12_fractures_test
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_fractures_test.csv 2021_12_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc


odir=2021_12_fractures_narrow
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_fractures_narrow.csv 2021_12_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc


odir=2021_12_all
grid=600
d=jib
n=24

python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 68:00:00 -n ${n} -g ${grid} -e ../uncertainty_quantification/ensemble_jib_all.csv 2021_12_init/state/jib_g600m_v1_RAGIS_id_INIT-TM-0.8-100-1.0_1980-1-1_1990-1-1.nc

# ISMIP 6 Runs

odir=2019_08_calib
grid=450
n=560

for d in gris; do  python historical.py --spatial_ts none -b no_bath -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell  -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_08_calib
n=224
grid=1000

for d in ismip6; do     python historical.py --spatial_ts none -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_08_calib
n=224
grid=1000

for d in ismip6; do     python historical.py --spatial_ts none -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2019_11_calib
n=224
grid=1000

for d in ismip6; do     python historical.py --spatial_ts ismip6 -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

# Hindcasts


odir=2020_09_hind
n=120
grid=600
for TCT in 250; do
    for d in jib; do python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${[d} --o_dir ${odir} --start 1985-1-1 --end 2000-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2020_09_sensitivity/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-${TCT}M-34S_1980-1-1_1985-1-1.nc ;
    done
    sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-${TCT}M-34S_1985-1-1_2000-1-1.sh
done






odir=2021_03_qaamerujup_h
n=48
grid=150

export PISM_PREFIX=~/pism-regional/bin; for d in qaamerujup; do     python historical.py -r 2 --stress_balance ssa+sia --spatial_ts basic --exstep monthly --dataset_version 4 -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2009-1-1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

export PISM_PREFIX=~/pism-regional/bin; for d in qaamerujup; do     python historical.py -r 2  --stress_balance ssa+sia --spatial_ts basic --exstep monthly --dataset_version 4 -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2009-1-1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv 2021_03_qaamerujup/state/qaamerujup_g150m_v4_id_CTRL_2008-1-1_2015-1-1.nc ; done

export PISM_PREFIX=~/pism-regional/bin; for d in qaamerujup; do     python historical.py -r 2  --stress_balance ssa+sia --spatial_ts basic --exstep 0:0.000001:0.0001 --dataset_version 4 -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2008-1-3 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv 2021_03_qaamerujup/state/qaamerujup_g150m_v4_id_CTRL_2008-1-1_2015-1-1.nc ; done


python historical.py --stress_balance blatter --spatial_ts none --dataset_version 3a -b wc -d gris --o_dir 2021_03_blatter --start 2008-1-1 --end 2015-1-1 -q debug -s chinook -w 1:00:00 -n 24 -g 9000 -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g9000m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc


python historical.py --options "'pc_type mg -stress_balance.blatter.Mz  65'" --stress_balance blatter --spatial_ts none --dataset_version 3a -b wc -d gris --o_dir 2021_03_blatter --start 2008-1-1 --end 2015-1-1 -q debug -s debug -w 1:00:00 -n 1 -g 18000 -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g9000m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc

odir=2021_03_init
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep daily --tsstep daily --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1985-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g450m_v1_RAGIS_id_0_0_50.nc ; done


odir=2021_03_melange
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep daily --tsstep daily --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1985-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2021_03_init/state/jib_g450m_v1_RAGIS_id_CTRL_1980-1-1_1985-1-1.nc  ; done

for file in ex_*; do
    ncks -O -L 2 -4 -d x,-208000.,-180000. -d y,-2285000.,-2255000. $file sm_$file;
    cdo -P 8 runmean,11 -fldmean sm_$file runmean_11days_fldmean_sm_$file;
done


odir=2021_04_theta
n=96
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1985-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_05_v1980_v3/state/gris_g900m_v1980v3_id_001_0_50.nc  ;
done


odir=2021_04_init
n=72
grid=900
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1988-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done

odir=2021_04_init
n=72
grid=900
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1988-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done

odir=2021_04_rm
n=72
grid=900
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b rm -d ${d} --o_dir ${odir} --start 1988-1-1 --end 1995-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done

odir=2021_04_rm
n=72
grid=900
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b rm -d ${d} --o_dir ${odir} --start 1988-1-1 --end 1995-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2021_04_rm/state/jib_g900m_v1_RAGIS_id_INIT-0.4-100-1.5_1980-1-1_1988-1-1.nc ;
done



odir=2021_04_rm
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b rm -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1988-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done


odir=2021_04_hind
n=72
grid=900
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1986-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2021_04_init/state/jib_g900m_v1_RAGIS_id_INIT-0.4-100-1.5_1980-1-1_1988-1-1.nc;
done

odir=2021_04_calib
n=72
grid=900
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b rm -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done


odir=2021_04_hc
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b rm -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2015-1-1 -q t2standard -s chinook -w 72:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done

odir=2021_04_calib
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b rm -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done


odir=2021_04_calib_wc
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b wc --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done


odir=2021_04_calib_rm
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done


odir=2021_04_init_rm
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done

odir=2021_04_init_rumpel
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done

odir=2021_04_calib_rumpel_0.75
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/calib_jib.csv 2021_04_init_rumpel/state/jib_g450m_v1_RAGIS_id_INIT-0.8-100-0.75_1980-1-1_1990-1-1.nc ;
done

odir=2021_04_calib_rumpel_1.00
n=120
grid=450
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/calib_jib.csv 2021_04_init_rumpel/state/jib_g450m_v1_RAGIS_id_INIT-0.8-100-1.00_1980-1-1_1990-1-1.nc ;
done


odir=2021_05_init
n=48
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2small -s chinook -w 36:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/init_jib.csv ../../pism-gris/calibration/2020_10_RAGIS/state/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc;
done

odir=2021_05_calib_1.00
n=48
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/calib_jib.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_1980-1-1_1990-1-1.nc  ;
done


odir=2021_05_calib_0.50
n=48
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/calib_jib.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-0.50_1980-1-1_1990-1-1.nc  ;
done

odir=2021_05_calib_1.00_blatter
n=48
grid=600
for d in jib; do     python historical.py --stress_balance blatter --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/calib_jib.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_1980-1-1_1990-1-1.nc  ;
done


odir=2021_02_blatter
n=96
grid=9000

for d in gris; do     python historical.py --stress_balance blatter --spatial_ts none --dataset_version 3a -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2009-1-1 -q t2standard -s chinook -w 1:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g9000m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2021_05_blatter
n=24
grid=9000

for d in gris; do     python historical.py --stress_balance blatter --spatial_ts none --dataset_version 3a -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2009-1-1 -q t2small -s chinook -w 1:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g9000m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2021_05_ocean
n=48
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_1980-1-1_1990-1-1.nc  ;
done

odir=2021_05_ocean_variability
n=48
grid=600
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1985-1-1 -q t2small -s chinook -w 4:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_1980-1-1_1990-1-1.nc  ;
    for id in CF; do
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/2021_05_ocean_variability/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_1985-1-1.sh
        python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2small -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv $odir/state/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_1985-1-1.nc;
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/2021_05_ocean_variability/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.sh
    done
done

odir=2021_06_ocean_variability
n=48
grid=600
a="1980-1-1"
e_1="1990-1-1"
e_2="2010-1-1"
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_1 -q t2small -s chinook -w 4:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a}_${e_1}.nc  ;
    for id in CTRL CF; do
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_1}.sh
        python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_2 -q t2small -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv $odir/state/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_1}.nc;
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_2}.sh
    done
done

odir=2021_07_ocean_variability
n=48
grid=600
a="1980-1-1"
e_1="2000-1-1"
e_2="2010-1-1"
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_1 -q t2small -s chinook -w 4:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a}_1990-1-1.nc  ;
    for id in CTRL CF {0..10}; do
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_1}.sh
        python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_2 -q t2small -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv $odir/state/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_1}.nc;
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_2}.sh
    done
done

odir=2021_07_calving_variability
n=48
grid=600
a="1980-1-1"
e_1="2000-1-1"
e_2="2010-1-1"
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a}_1990-1-1.nc  ;
#    for id in VCM-0.4-1.00 VCM-0.4-1.25 VCM-0.4-1.50 VCM-0.5-1.00 VCM-0.5-1.25 VCM-0.5-1.50 VCM-0.6-1.00 VCM-0.6-1.25 VCM-0.6-1.50 VCM-0.7-1.00  VCM-0.8-1.00 ; do
    for id in  VCM-0.8-1.00-200 ; do
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_CALV-${id}_${a}_${e_1}.sh
        python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_2 -q t2small -s chinook -w 16:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving.csv ${odir}/state/jib_g600m_v1_RAGIS_id_CALV-${id}_${a}_${e_1}.nc;
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_CALV-${id}_${a}_${e_2}.sh
    done
done


odir=2021_07_calving_variability
n=48
grid=600
a="1980-1-1"
e_1="2000-1-1"
e_2="2010-1-1"
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a}_1990-1-1.nc  ;
#    for id in VCM-0.4-1.00 VCM-0.4-1.25 VCM-0.4-1.50 VCM-0.5-1.00 VCM-0.5-1.25 VCM-0.5-1.50 VCM-0.6-1.00 VCM-0.6-1.25 VCM-0.6-1.50 VCM-0.7-1.00  VCM-0.8-1.00 ; do
    for id in  VCM-0.7-1.00  VCM-0.8-1.00 ; do
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_CALV-${id}_${a}_${e_1}.sh
        python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_2 -q t2small -s chinook -w 16:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving.csv ${odir}/state/jib_g600m_v1_RAGIS_id_CALV-${id}_${a}_${e_1}.nc;
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_CALV-${id}_${a}_${e_2}.sh
    done
done

odir=2021_07_calving_calib
n=48
grid=600
a="1980-1-1"
e_1="2000-1-1"
e_2="1985-8-1"
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_calib.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a}_1990-1-1.nc  ;
    for id in CALV-0.5-1.00-200 CALV-0.5-1.25-200  CALV-0.5-1.50-200  CALV-0.5-1.00-250 CALV-0.5-1.00-250 CALV-0.5-1.25-250 CALV-0.5-1.50-250; do
#     for id in CALV-0.6-1.00-200 CALV-0.7-1.00-200 CALV-0.8-1.00-200 CALV-0.6-1.25-200 CALV-0.7-1.25-200 CALV-0.8-1.25-200 CALV-0.6-1.50-200 CALV-0.7-1.50-200 CALV-0.8-1.50-200 CALV-0.6-1.00-250 CALV-0.7-1.00-250 CALV-0.8-1.00-250 CALV-0.6-1.25-250 CALV-0.7-1.25-250 CALV-0.8-1.25-250 CALV-0.6-1.50-250 CALV-0.7-1.50-250 CALV-0.8-1.50-250; do
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_${id}_${a}_${e_1}.sh
        # python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_2 -q t2small -s chinook -w 5:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_calib.csv ${odir}/state/jib_g600m_v1_RAGIS_id_CALV-${id}_${a}_${e_1}.nc;
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_${id}_${a}_${e_2}.sh
    done
done

odir=2021_07_calving_calib
n=48
grid=600
a="1980-1-1"
e_1="2000-1-1"
e_2="1985-8-1"
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_1 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_calib.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a}_1990-1-1.nc  ;
    for id in CALV-0.5-1.00-200 CALV-0.5-1.25-200  CALV-0.5-1.50-200  CALV-0.5-1.00-250 CALV-0.5-1.25-250 CALV-0.5-1.50-250; do
#     for id in CALV-0.6-1.00-200 CALV-0.7-1.00-200 CALV-0.8-1.00-200 CALV-0.6-1.25-200 CALV-0.7-1.25-200 CALV-0.8-1.25-200 CALV-0.6-1.50-200 CALV-0.7-1.50-200 CALV-0.8-1.50-200 CALV-0.6-1.00-250 CALV-0.7-1.00-250 CALV-0.8-1.00-250 CALV-0.6-1.25-250 CALV-0.7-1.25-250 CALV-0.8-1.25-250 CALV-0.6-1.50-250 CALV-0.7-1.50-250 CALV-0.8-1.50-250; do
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_${id}_${a}_${e_1}.sh
        python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_2 -q t2small -s chinook -w 5:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_calib.csv ${odir}/state/jib_g600m_v1_RAGIS_id_${id}_${a}_${e_1}.nc;
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_${id}_${a}_${e_2}.sh
    done
done


odir=2021_08_calving_threshold
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1985-1-1"
e_1="1985-1-1"
e_2="2010-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_threshold.csv 2021_07_calving_calib/state/jib_g600m_v1_RAGIS_id_CALV-0.5-1.25-250_1980-1-1_1985-8-1.nc;

odir=2021_08_healing
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1985-1-1"
e_1="1985-1-1"
e_2="2010-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_threshold.csv 2021_07_calving_calib/state/jib_g600m_v1_RAGIS_id_CALV-0.5-1.25-250_1980-1-1_1985-8-1.nc;

odir=2021_08_healing_var
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1985-1-1"
e_1="1985-1-1"
e_2="2010-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_threshold.csv 2021_07_calving_calib/state/jib_g600m_v1_RAGIS_id_CALV-0.5-1.25-250_1980-1-1_1985-8-1.nc;

odir=2021_08_fractures
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1980-1-1"
e_1="1990-1-1"
e_2="2010-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_fractures.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a_1}_${e_1}.nc;

odir=2021_09_fractures
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1980-1-1"
e_1="1990-1-1"
e_2="2010-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 28:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_fractures.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a_1}_${e_1}.nc;


odir=2021_09_uq_10
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1980-1-1"
e_1="1990-1-1"
e_2="1986-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 8:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_uq_samples_10.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a_1}_${e_1}.nc;


odir=2021_09_uq_borstadt_2
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1980-1-1"
e_1="1990-1-1"
e_2="1988-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ensemble_jib_borstadt_2_uq_samples_10.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a_1}_${e_1}.nc;

odir=2021_09_uq_borstadt_20
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1980-1-1"
e_1="1990-1-1"
e_2="1988-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 16:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ensemble_jib_borstadt_samples_20.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a_1}_${e_1}.nc;


odir=2021_10_variability
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1980-1-1"
e_1="1990-1-1"
e_2="1988-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 16:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ensemble_jib_variability_16.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a_1}_${e_1}.nc;



odir=2021_10_hind
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1980-1-1"
e_1="1990-1-1"
e_2="2010-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ensemble_jib_borstadt_samples_20.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_${a_1}_${e_1}.nc;

odir=2021_08_calving_threshold_blatter
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1985-1-1"
e_1="1985-1-1"
e_2="2010-1-1"
python historical.py --stress_balance blatter --hydrology routing --spatial_ts standard --exstep daily --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q t2small -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_threshold.csv 2021_07_calving_calib/state/jib_g600m_v1_RAGIS_id_CALV-0.5-1.25-250_1980-1-1_1985-8-1.nc;


odir=2021_07_calving_threshold
n=24
grid=600
d=jib
a_1="1980-1-1"
a_2="1985-1-1"
e_1="1985-1-1"
e_2="2010-1-1"
python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a_2 --end $e_2 -q normal -s pleiades_broadwell -w 00:25:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_calving_threshold.csv 2021_07_calving_calib/state/jib_g600m_v1_RAGIS_id_CALV-0.5-1.25-250_1980-1-1_1985-8-1.nc;


odir=2021_06_test
n=48
grid=600
a="1960-1-1"
e_1="1990-1-1"
e_2="2010-1-1"
for d in jib; do
    python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_1 -q t2small -s chinook -w 4:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv 2021_05_init/state/jib_g600m_v1_RAGIS_id_INIT-0.8-100-1.00_1980-1-1_${e_1}.nc  ;
    for id in CTRL CF; do
        sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_1}.sh
        python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d ${d} --o_dir ${odir} --start $a --end $e_2 -q t2small -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/jib_ocean_variability.csv $odir/state/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_1}.nc;
        # sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_${a}_${e_2}.sh
    done
done



odir=2021_04_calib_1985
mkdir -p $odir/glaciers/scalar
cd $odir/spatial
for file in ex_jib*.nc; do
 python ~/base/gris-analysis/basins/extract_glacier.py --ugid 225 --epsg 3413 --o_dir ../glaciers --shape_file ~/base/gris-analysis/basins/Greenland_Basins_PS_v1.4.2_1980.shp ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-0_1980-1-1_2010-1-1.nc  $file
done
cd ../glaciers
for id in INIT1 INIT2 INIT3 INIT4; do
cdo -L -O expr,"subshelf_melt_rate=basal_mass_flux_floating/-1000.0" -fldmean -ifthen ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_${id}_1985-1-1/ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_${id}_1985-1-1.nc ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_${id}_1985-1-1/ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_${id}_1985-1-1.nc scalar/subshelf_melt_rate_ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_${id}_1985-1-1.nc    
done


cdo -L -O expr,"subshelf_melt_rate=basal_mass_flux_floating/-1000.0" -fldmean -ifthen ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_INIT1_1985-1-1/ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_INIT1_1985-1-1.nc ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_INIT1_1985-1-1/ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_INIT1_1985-1-1.nc scalar/subshelf_melt_rate_ugid_0_JIB_tongue_ex_jib_g600m_v1_RAGIS_id_INIT1_1985-1-1.nc


odir=2021_05_ocean_variability
mkdir -p $odir/glaciers/scalar
cd $odir/spatial
for file in ex_jib*2010*.nc; do
 python ~/base/gris-analysis/basins/extract_glacier.py --ugid 225 --epsg 3413 --o_dir ../glaciers --shape_file ~/base/gris-analysis/basins/Greenland_Basins_PS_v1.4.2_1980.shp $file
done
# cd ../glaciers
for id in {0..9}; do
cdo -L -O setattribute,subshelf_melt_rate@units="m yr-1" -aexpr,"subshelf_melt_rate=basal_mass_flux_floating/-1000.0" -fldmean -ifthen ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1/ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1/ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc scalar/fldmean_ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc
done
for id in {0..9}; do
cdo -L -O setattribute,subshelf_melt_rate@units="m yr-1" -aexpr,"subshelf_melt_rate=basal_mass_flux_floating/-1000.0" -fldsum -ifthen ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1/ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1/ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc scalar/fldsum_ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc
done


for id in {0..9} CTRL; do
extract_profiles.py -v velsurf_mag --srs epsg:3413 ~/Google\ Drive/My\ Drive/Projects/jib-breakup/data/shape_files/joughin-gps-points.shp ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1/ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc points_ugid_225_Jakobshavn_Isbrae_ex_jib_g600m_v1_RAGIS_id_OCEAN-VAR-${id}_1980-1-1_2010-1-1.nc
done

odir=2021_12_blatter
n=24
grid=4500

for d in gris; do     python historical.py --stress_balance blatter --spatial_ts standard --exstep daily --dataset_version 3a -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2009-1-1 -q t2small -s chinook -w 1:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2021_08_blatter
n=72
grid=1800

for d in gris; do     python historical.py --stress_balance blatter --spatial_ts standard --exstep daily --dataset_version 3a -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2008-6-30 -q t2standard -s chinook -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2021_08_blatter_adaptive_10
n=72
grid=1800

for d in gris; do     python historical.py --stress_balance blatter --spatial_ts standard --exstep daily --dataset_version 3a -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2008-6-30 -q t2standard -s chinook -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2021_08_blatter_adaptive_25
n=72
grid=1800

for d in gris; do     python historical.py --stress_balance blatter --spatial_ts standard --exstep daily --dataset_version 3a -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2008-6-30 -q t2standard -s chinook -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2021_08_blatter_adaptive_25
n=24
grid=4500

for d in gris; do     python historical.py --stress_balance blatter --spatial_ts standard --exstep daily --dataset_version 3a -b wc -d ${d} --o_dir ${odir} --start 2008-1-1 --end 2008-6-30 -q t2standard -s chinook -w 48:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/ctrl.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g${grid}m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done
