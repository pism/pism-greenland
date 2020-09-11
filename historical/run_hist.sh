

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

odir=2019_07_hind
n=224
grid=1000
for d in ismip6; do     python historical.py --spatial_ts basic --exstep monthly -b wc -d ${d} --o_dir ${odir} --end 2015-1-1 -q long -s pleiades_broadwell -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_discharge_given.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2020_05_gimp
n=240
grid=900
for d in gris; do     python historical.py --spatial_ts basic --exstep monthly --tsstep daily -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1998-1-1 -q t2standard -s chinook -w 60:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2020_05_gimp
n=240
grid=900
for d in gris; do     python historical.py --spatial_ts basic --exstep monthly --tsstep daily -b wc -d ${d} --o_dir ${odir} --start 1998-1-1 --end 2008-1-1 -q t2standard -s chinook -w 60:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv /import/c1/ICESHEET/ICESHEET/crios2pism/historical/2020_05_gimp/state/gris_g900m_v3a_id_HH5-CNRM-CM6-MM-VCM0.70-100M-34S_1980-1-1_1998-1-1.nc; done

# RELAX
grid=600
for d in jib; do     python historical.py --spatial_ts basic --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1981-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc ; done

for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1986-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv 2020_05_gimp/state/jib_g600m_v1980_id_MAR-CNRM-CM6-MM-VCM0.70-100M-34S_1980-1-1_1981-1-1.nc ; done

odir=2020_05_reg
n=240
grid=900
for d in gris; do     python historical.py --spatial_ts basic --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1982-1-1 --end 2008-1-1 -q t2standard -s chinook -w 60:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2020_05_v1980/state/gris_g900m_v1980_id_001_0_50.nc; done

odir=2020_05_reg_v3
n=240
grid=900
for d in gris; do     python historical.py --spatial_ts basic --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1982-1-1 --end 2008-1-1 -q t2standard -s chinook -w 60:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2020_05_v1980_v3/state/gris_g900m_v1980_id_001_0_50.nc ; done


# RELAX
grid=600
for d in jib; do     python historical.py --spatial_ts basic --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1981-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2020_05_v1980_v3/state/gris_g900m_v1980v3_id_001_0_50.nc ; done

for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1986-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv 2020_05_reg_v3/state/jib_g600m_v1980_id_MAR-CNRM-CM6-MM-VCM0.70-100M-34S_1980-1-1_1981-1-1.nc ; done

for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1986-1-1 --end 2008-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv 2020_05_reg_v3/state/jib_g600m_v1980_id_MAR-CNRM-CM6-MM-VCM0.60-175M-34S_1980-1-1_1986-1-1.nc ; done


# RELAX
odir=2020_09_relax
n=120
grid=600
for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2020_05_v1980_v3/state/gris_g900m_v1980v3_id_001_0_50.nc ; done


# RELAX
odir=2020_09_relax_routing
n=120
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_05_v1980_v3/state/gris_g900m_v1980v3_id_001_0_50.nc ; done


# RELAX
odir=2020_09_relax_p
n=120
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1990-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2020_09_relax_p/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-100M-34S_1980-1-1_1990-1-1.nc ; done

for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1985-1-1 --end 2000-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2020_09_relax_p/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-250M-34S_1980-1-1_1990-1-1.nc  ; done


odir=2020_09_init
n=120
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1985-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv ../../pism-gris/calibration/2020_05_v1980_v3/state/gris_g900m_v1980v3_id_001_0_50.nc ; done

odir=2020_09_sensitivity
n=120
grid=600
for d in jib; do     python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1985-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2020_09_init/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-100M-34S_1980-1-1_1985-1-1.nc ; done

odir=2020_09_hind
n=120
grid=600
for TCT in 250; do
    for d in jib; do python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1985-1-1 --end 2000-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2020_09_sensitivity/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-${TCT}M-34S_1980-1-1_1985-1-1.nc ;
    done
    sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-${TCT}M-34S_1985-1-1_2000-1-1.sh
done


odir=2020_09_hind
n=120
grid=600
for TCT in 250; do
    for d in jib; do python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2020_09_sensitivity/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-${TCT}M-34S_1980-1-1_1985-1-1.nc ;
    done
    sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-${TCT}M-34S_1980-1-1_2010-1-1.sh
done

odir=2020_09_hind
n=120
grid=600
for TCT in 325; do
    for d in jib; do python historical.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep monthly --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2010-1-1 -q t2standard -s chinook -w 12:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_jib.csv 2020_09_sensitivity/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.65-250M-34S_1980-1-1_1985-1-1.nc ;
    done
    sbatch /import/c1/ICESHEET/ICESHEET/crios2pism/historical/${odir}/run_scripts/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.50-${TCT}M-34S_1980-1-1_2010-1-1.sh
done

for file in ex_jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0*2010-1-1.nc; do
    echo $file
    ncks -4 -d x,-208000.,-180000. -d y,-2285000.,-2255000. $file sm_$file
    cdo -O -L -f nc4 -P 8 runmean,13 sm_$file runmean_13mo_sm_$file
    cdo -O -L runmean,13 -fldmean -ifthen sm_$file sm_$file  runmean_13mo_fldmean_sm_$file
done


odir=2020_05_gimp
n=240
grid=1000
for d in ismip6; do     python historical.py --spatial_ts basic --exstep monthly --tsstep daily -b wc -d ${d} --o_dir ${odir} --start 1982-1-1 --end 2008-1-1 -q t2standard -s chinook -w 60:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc; done

odir=2020_05_tct
n=120
grid=600
for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1982-1-1 --end 1988-1-1 -q t2standard -s chinook -w 10:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2020_05_v1980/state/gris_g900m_v1980_id_001_0_50.nc; done

odir=2020_05_tct
n=120
grid=600
for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1988-1-1 --end 2008-1-1 -q t2standard -s chinook -w 10:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv 2020_05_tct/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.80-50M_1982-1-1_1988-1-1.nc; done

odir=2020_05_hind_v2
n=120
grid=600
for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1986-1-1 -q t2standard -s chinook -w 10:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2020_05_v1980_v2/state/gris_g900m_v1980_id_001_0_50.nc; done

odir=2020_05_reg_v2
n=360
grid=900
for d in gris; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1986-1-1 -q t2standard -s chinook -w 72:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2020_05_v1980_v2/state/gris_g900m_v1980_id_001_0_50.nc; done


odir=2020_05_hind_c
n=120
grid=600
for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1986-1-1 --end 1998-1-1 -q t2standard -s chinook -w 10:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv 2020_05_hind_c/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.70-100M-34S_1980-1-1_1986-1-1.nc  ; done

odir=2020_05_hind_c
n=120
grid=600
for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1998-1-1 --end 2008-1-1 -q t2standard -s chinook -w 10:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv 2020_05_hind_c/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.70-100M-34S_1986-1-1_1998-1-1.nc  ; done


odir=2020_05_hind_v2
n=120
grid=600

for d in jib; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1998-1-1 --end 2008-1-1 -q t2standard -s chinook -w 10:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv  2020_05_hind/state/jib_g600m_v1980_id_HH5-CNRM-CM6-MM-VCM0.80-50M-34S_1982-1-1_1998-1-1.nc ; done

odir=2020_05_hind_c
n=240
grid=900
for d in gris; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 1998-1-1 -q t2standard -s chinook -w 60:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv /import/c1/ICESHEET/ICESHEET/crios2pism/historical/2020_05_reg_v2/state/gris_g900m_v1980_id_HH5-CNRM-CM6-MM-VCM0.70-100M-34S_1980-1-1_1986-1-1.nc  ; done

odir=2020_05_hind_c
n=240
grid=900

for d in gris; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily --dataset_version 1980 -b wc -d ${d} --o_dir ${odir} --start 1998-1-1 --end 2008-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv  /import/c1/ICESHEET/ICESHEET/crios2pism/historical/2020_05_hind_c/state/gris_g900m_v1980_id_HH5-CNRM-CM6-MM-VCM0.70-100M-34S_1980-1-1_1998-1-1.nc ; done



/import/c1/ICESHEET/ICESHEET/crios2pism/historical/2020_05_reg_v2/state/gris_g900m_v1980_id_HH5-CNRM-CM6-MM-VCM0.70-100M-34S_1980-1-1_1986-1-1.nc

odir=2020_08_fd
n=240
grid=900

for d in gris; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily  -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2015-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc ; done

odir=2020_08_fem
n=240
grid=900

for d in gris; do     python historical.py --spatial_ts standard --exstep monthly --tsstep daily  -b wc -d ${d} --o_dir ${odir} --start 1980-1-1 --end 2015-1-1 -q t2standard -s chinook -w 24:00:00 -n ${n} -g ${grid} -e ../uncertainty_qunatification/historical_hirham.csv ../../pism-gris/calibration/2017_06_vc/state/gris_g900m_flux_v3a_no_bath_sia_e_1.25_sia_n_3_ssa_n_3.25_ppq_0.6_tefo_0.02_calving_vonmises_calving_0_100.nc ; done

