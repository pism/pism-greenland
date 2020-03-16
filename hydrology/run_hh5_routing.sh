#!/bin/sh
#SBATCH --partition=t2small
#SBATCH --ntasks=48
#SBATCH --tasks-per-node=24
#SBATCH --time=16:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=pism.%j

module list

cd $SLURM_SUBMIT_DIR

# Generate a list of compute node hostnames reserved for this job,
# this ./nodes file is necessary for slurm to spawn mpi processes
# across multiple compute nodes
srun -l /bin/hostname | sort -n | awk '{print $2}' > ./nodes_$SLURM_JOBID

ulimit -l unlimited
ulimit -s unlimited
ulimit


# stop if a variable is not defined
set -u
# stop on errors
set -e

# path to the input directory (input data sets are contained in this directory)
input_dir="/import/c1/ICESHEET/ICESHEET/crios2pism"
# output directory
output_dir="//import/c1/ICESHEET/ICESHEET/crios2pism/hydrology/2020_03_routing"

# create required output directories
for each in $output_dir;
do
  mkdir -p $each
done

GRID=$1

mpiexec -n 48 $HOME/pism/bin/pismr \
        -i ../data_sets/bed_dem/pism_Greenland_${GRID}m_mcb_jpl_v4_ctrl.nc \
        -o_size none \
        -bootstrap \
        -Mz 3 \
        -time_file ../data_sets/runoff/DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_4500m_MM.nc \
        -hydrology routing \
        -hydrology.tillwat_max 0 \
        -stress_balance none \
        -energy none \
        -hydrology.surface_input.file ../data_sets/runoff/DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_${GRID}m_MM.nc \
        -extra_times monthly \
        -extra_vars bwat,tillwat,hydrology_fluxes,subglacial_water_input_rate,subglacial_water_flux_mag \
        -extra_file $output_dir/ex_g${GRID}m_water_routing_DMI-HIRHAM5_GL2_ERAI_1980_2016_MM.nc \
         > $output_dir/job.$SLURM_JOBID 2>&1

ncks -O -4 -L 2 $output_dir/ex_g${GRID}m_water_routing_DMI-HIRHAM5_GL2_ERAI_1980_2016_MM.nc $output_dir/ex_g${GRID}m_water_routing_DMI-HIRHAM5_GL2_ERAI_1980_2016_MM.nc
