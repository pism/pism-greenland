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

ofile_dm=DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_${grid}M_DM.nc
ofile_mm=DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_EPSG3413_${grid}M_MM.nc
tmp_file=tmp_${ofile_dm}
#
# "DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_DM.nc" in mm / day with rho = 1000 kg / m3 we get kg / m2 / day
# "DMI-HIRHAM5_GL2_ERAI_2015_2016_MRROS_DM.nc" claims to be in kg / m2 / s but this doesn't make sense,
# it's in kg m-2 day-1
cdo -L -O -f nc4 --reduce_dim sellevel,0 -mergetime -setmisstoc,0 -setattribute,mrros@units="kg m-2 day-1" ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_DM.nc -setattribute,mrros@units="kg m-2 day-1" ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_2015_2016_MRROS_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_DM.nc
adjust_timeline.py -p daily -a 1980-1-1 -d 1980-1-1 DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_DM.nc
create_greenland_epsg3413_grid.py -g ${grid} $grid_file
nc2cdo.py $grid_file
cdo remapbil,${grid_file} -setgrid,../../util/rotated_grid.txt DMI-HIRHAM5_GL2_ERAI_1980_2016_MRROS_DM.nc $tmp_file
cdo -O -L setmisstoc,0 -setattribute,water_input_rate@units="kg m-2 day-1" -chname,mrros,water_input_rate $tmp_file $ofile_dm
cdo monmean $ofile_dm $ofile_mm
rm $tmp_file
