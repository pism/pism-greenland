#!/bin/bash
#SBATCH --partition=analysis
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --output=pism.%j
#SBATCH --mem=214G

cd $SLURM_SUBMIT_DIR

odir=$1
file=$2

mkdir -p ${odir}/io

# Spatial files
cd  ${odir}/spatial/

cdo -f nc4 -z zip_2 ifnotthen ../../../data_sets/basin_masks/ugid_225_Jakobshavn_Isbrae_mask_epsg3413_g600m.nc ex_$file masked_ex_$file
cdo fldsum  masked_$file fldsum_masked_$file
cdo fldmean masked_$file mean_masked_$file

extract_interface.py -t ice_ocean -o ../io/io_masked_ex_$file masked_ex_$file

# State files
cd ../state

cdo -f nc4 -z zip_2 ifnotthen ../../../data_sets/basin_masks/ugid_225_Jakobshavn_Isbrae_mask_epsg3413_g600m.nc $file masked_$file
extract_interface.py -t ice_ocean -o ../io/io_masked_$file masked_$file

