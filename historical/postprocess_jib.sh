#!/bin/bash
#SBATCH --partition=analysis
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --output=pism.%j
#SBATCH --mem=214G

cd $SLURM_SUBMIT_DIR

odir=$1
file=$2
grid=$3

mkdir -p ${odir}/io
mkdir -p ${odir}/profiles

# Spatial files
cd  ${odir}/spatial/

cdo -f nc4 -z zip_2 aexpr,"total_grounding_line_flux=grounding_line_flux*${grid}^2/1e12;" -ifnotthen ../../../data_sets/basin_masks/ugid_225_Jakobshavn_Isbrae_mask_epsg3413_g${grid}m.nc ex_${file}.nc masked_ex_${file}.nc
ncatted -a units,total_grounding_line_flux,o,c,"Gt year-1" masked_ex_${file}.nc

for var in velsurf_mag thk; do
cdo -f nc4 -z zip_2 seldate,1985-7-16 masked_ex_${file}.nc ${var}_masked_ex_${file}_1985-7-16.nc
cdo -f nc4 -z zip_2 seldate,2010-7-16 masked_ex_${file}.nc ${var}_masked_ex_${file}_1985-7-16.nc
done

cdo fldsum  masked_ex_${file}.nc fldsum_masked_ex_${file}.nc
cdo fldmean masked_ex_${file}.nc fldmean_masked_ex_${file}.nc

extract_interface.py -t ice_ocean -o ../io/io_masked_ex_${file}.nc masked_ex_${file}.nc
extract_profiles.py -v velsurf_mag --srs epsg:3413 /import/c1/ICESHEET/ICESHEET/crios2pism/data_sets/shape_files/joughin-gps-points.shp ex_${file}.nc ../profiles/gps_stations_ex_${file}.nc

# State files
cd ../state

cdo -f nc4 -z zip_2 ifnotthen ../../../data_sets/basin_masks/ugid_225_Jakobshavn_Isbrae_mask_epsg3413_g${grid}m.nc ${file}.nc masked_${file}.nc
extract_interface.py -t ice_ocean -o ../io/io_masked_${file}.gkpg masked_${file}.nc

