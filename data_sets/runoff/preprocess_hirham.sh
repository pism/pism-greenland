#!/bin/bash

prefix="http://prudence.dmi.dk/data/temp/RUM/HIRHAM/GREENLAND/ERAI/DAILY"

gz_dir=gz_files
extracted_dir=extracted_files
grid=4500
grid_file=g${grid}m.nc

for dir in "$gz_dir"; do
    mkdir -p $dir
done

for file in "DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_DM.nc" "DMI-HIRHAM5_GL2_ERAI_2015_2016_MRROS_DM.nc"; do
    wget -nc ${prefix}/${file}.gz -P ${gz_dir}/
        gunzip -c ${gz_dir}/${file}.gz > ${extracted_dir}/${file}
done

ofile_dm=DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_DM_EPSG3413_${grid}m.nc
ofile_ydm=DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_YDM_EPSG3413_${grid}m.nc
tmp_file=tmp_${ofile_dm}
# The two files have different input units, so we ignore the newer one for now.
#cdo -O mergetime ${gz_dir}/*.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_DM.nc
cp ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_DM.nc
create_greenland_epsg3413_grid.py -g 4500 $grid_file
cdo remapbil,${grid_file} -setgrid,../../util/rotated_grid.txt DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_DM.nc $tmp_file
cdo -O setmisstoc,0 -chname,mrros,water_input_rate $tmp_file $ofile_dm
adjust_timeline.py -p daily -a 1980-1-1 -d 1980-1-1 $ofile_dm
cdo ydaymean $ofile_dm $ofile_ydm
rm $tmp_file
