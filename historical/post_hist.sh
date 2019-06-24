#!/bin/bash
#SBATCH --partition=analysis
#SBATCH --ntasks=7
#SBATCH --tasks-per-node=7
#SBATCH --time=48:00:00
#SBATCH --output=pism.%j
#SBATCH --mem=214G

cd $SLURM_SUBMIT_DIR

odir=$1
exp=$2
grid=$3
domain=$4

mkdir -p ${odir}/io

extract_interface.py -t ice_ocean -o ${odir}/io/io_${domain}_g${grid}m_v3a_id_${exp}_2008-1-1_2015-1-1.nc ${odir}/state/${domain}_g${grid}m_v3a_id_${exp}_2008-1-1_2015-1-1.nc 
ncks -O -4 -L 3 ${odir}_tmp/ex_${domain}_g${grid}m_v3a_id_${exp}_2008-1-1_2015-1-1.nc ${odir}/spatial/ex_${domain}_g${grid}m_v3a_id_${exp}_2008-1-1_2015-1-1.nc
ncks -O -4 -L 3 ${odir}/state/${domain}_g${grid}m_v3a_id_${exp}_2008-1-1_2015-1-1.nc  ${odir}/state/${domain}_g${grid}m_v3a_id_${exp}_2008-1-1_2015-1-1.nc 
