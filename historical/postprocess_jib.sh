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
ensfile=$4

mkdir -p ${odir}/io
mkdir -p ${odir}/profiles
mkdir -p ${odir}/processed

# Spatial files
cd  ${odir}/spatial/

cdo -f nc4 -z zip_2 setctomiss,0 -aexpr,"total_grounding_line_flux=grounding_line_flux*${grid}^2/1e12;" -ifnotthen ../../../data_sets/basin_masks/ugid_225_Jakobshavn_Isbrae_mask_epsg3413_g${grid}m.nc ex_${file}.nc ../processed/masked_ex_${file}.nc
ncatted -a units,total_grounding_line_flux,o,c,"Gt year-1" ../processed/masked_ex_${file}.nc

# for var in mask velsurf_mag thk; do
#     cdo -f nc4 -z zip_2 seldate,1985-7-16 -selvar,$var ../processed/masked_ex_${file}.nc ../processed/${var}_masked_ex_${file}_1985-7-16.nc
#     cdo -f nc4 -z zip_2 seldate,2009-7-16 -selvar,$var ../processed/masked_ex_${file}.nc ../processed/${var}_masked_ex_${file}_2009-7-16.nc
# done

cdo fldsum ../processed/masked_ex_${file}.nc ../processed/fldsum_masked_ex_${file}.nc
cdo fldmean ../processed/masked_ex_${file}.nc ../processed/fldmean_masked_ex_${file}.nc

# extract_interface.py --ensemble_file ../../../uncertainty_quantification/$ensfile -e 3413 -t ice_ocean -o ../io/io_masked_ex_${file}_1985-7-16.gpkg ../processed/mask_masked_ex_${file}_1985-7-16.nc

# extract_profiles.py -v velsurf_mag --srs epsg:3413 /import/c1/ICESHEET/ICESHEET/crios2pism/data_sets/shape_files/joughin-gps-points.shp ex_${file}.nc ../profiles/gps_stations_ex_${file}.nc

# ~/base/pypismtools/scripts/extract_profiles.py --srs epsg:3413 ~/base/gris-analysis/flux-gates/jakobshavn-flowline-50m.shp ex_${file}.nc ../profiles/flowline_ex_${file}.nc

# cd ../state/
# extract_interface.py --ensemble_file ../../../uncertainty_quantification/$ensfile -e 3413 -t ice_ocean -o ../io/io_${file}_2010-1-1.gpkg ${file}.nc

# cd ../../
# python ../util/create_flowline_animation.py -b ${odir} ${odir}/profiles/flowline_ex_${file}.nc
