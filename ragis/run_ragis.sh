n=68
grid=1800
odir=2022_08_init

PISM_PREFIX=$WORK/local/pism-dev/bin/ python3 hindcast.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d gris --o_dir $SCRATCH/${odir} --start 1980-1-1 --end 1990-1-1 -q normal -s stampede2 -w 4:00:00 -n 68 -g ${grid} -e $WORK/crios2pism/uncertainty_quantification/ensemble_gris_ragis_calving_lhs_20_w_posterior.csv $WORK/crios2pism/data_sets/initial_states/gris_g${grid}m_v1_RAGIS_id_0_0_50.nc

n=68
grid=1800
odir=2022_08_calving

for id in {0..19}; do
PISM_PREFIX=$WORK/local/pism-dev/bin/ python3 hindcast.py --hydrology routing --spatial_ts standard --exstep monthly --tsstep daily -b rm --dataset_version 1_RAGIS -d gris --o_dir $SCRATCH/${odir} --start 1980-1-1 --end 2020-1-1 -q normal -s stampede2 -w 8:00:00 -n 68 -g ${grid} -e $WORK/crios2pism/uncertainty_quantification/ensemble_gris_ragis_calving_lhs_20_w_posterior.csv $SCRATCH/2022_08_init/state/gris_g1800m_v1_RAGIS_id_${id}_1980-1-1_1990-1-1.nc

sbatch $SCRATCH/${odir}/run_scripts/gris_g1800m_v1_RAGIS_id_${id}_1980-1-1_2020-1-1.sh
done