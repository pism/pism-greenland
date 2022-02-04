#!/bin/bash

hh5_dir=hirham5_zip
hh5_nc_dir=hirham5_nc
mkdir -p $hh5_dir
mkdir -p $hh5_nc_dir

for year in {1980..2021}; do
    wget -nc -P $hh5_dir http://ensemblesrt3.dmi.dk/data/prudence/temp/nichan/Daily2D_GrIS/${year}.zip
    unzip -u ${hh5_dir}/${year}.zip
    cdo -O -f nc4 -z zip_2 -r settbounds,day -setattribute,gld@units="kg m-2 day-1",rogl@units="kg m-2 day-1",snfall@units="kg m-2 day-1",rainfall@units="kg m-2 day-1" -selvar,gld,tas,rogl,snfall,rainfall -mergetime ${year}/Daily2D_*.nc ${hh5_nc_dir}/${year}.nc
done


start_year=2019
end_year=2020
for year in 2019; do
    cdo -O -f nc4 -r settbounds,day -mergetime Daily2D_GrIS/${year}/Daily*.nc DMI-HIRHAM5_${year}.nc
done
cdo -O -f nc4 mergetime DMI-HIRHAM5_${year}.nc DMI-HIRHAM5_${start_year}_${end_year}.nc

cdo -O -P 6 -f nc4   setmisstoc,0 -selvar,climatic_mass_balance,ice_surface_temp,precipitation,water_input_rate -setattribute,climatic_mass_balance@units="kg m-2 day-1" -chname,rogl,water_input_rate -chname,gld,climatic_mass_balance -chname,tas,ice_surface_temp -chname,rainfall,precipitation -remapbil,../../grids/gris_ext_g4500m.txt -setgrid,../../grids/rotated_grid.txt  DMI-HIRHAM5_${year}.nc test.nc


 
# Download precipitation and 2-m air temperature
for var in PR TAS; do
    # for file in DMI-HIRHAM5_GL2_ERAI_1980_2014_${var}_DM.nc DMI-HIRHAM5_GL2_ERAI_2015_2016_${var}_DM.nc; do 
    #     wget -nc ${prefix}/${file}.gz -P ${gz_dir}/
    #     gunzip -c ${gz_dir}/${file}.gz > ${extracted_dir}/${file}
    # done
    # remove vertical levels and merge files
    # adjust_timeline.py -p daily -a 1980-1-1 -d 1980-1-1 ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2014_${var}_DM.nc
    # adjust_timeline.py -p daily -a 2015-1-1 -d 1980-1-1 ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_2015_2016_${var}_DM.nc    
    cdo -O -f nc4 --reduce_dim mergetime ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2014_${var}_DM.nc ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_2015_2016_${var}_DM.nc ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2016_${var}_DM.nc
done
exit
# # Download climatic mass balance
file=DMI-HIRHAM5_GL2_ERAI_1980_2016_gld_DM.nc
wget -nc http://prudence.dmi.dk/data/temp/RUM/HIRHAM/GREENLAND/ERAI/DAILY/${file}.gz -P ${gz_dir}/
gunzip -c ${gz_dir}/${file}.gz > ${extracted_dir}/${file}

# Merge precipitation, 2-m air temperature, and climatic mass balance
cdo -O -f nc4 -z zip_2 merge  ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2016_gld_DM.nc ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2016_TAS_DM.nc ${extracted_dir}/DMI-HIRHAM5_GL2_ERAI_1980_2016_PR_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_DM.nc

# Create target grid
create_greenland_epsg3413_grid.py -g ${grid} ${grid_file}
nc2cdo.py ${grid_file}

# Adjust the time line
adjust_timeline.py -p daily -a 1980-1-1 -d 1980-1-1 DMI-HIRHAM5_GL2_ERAI_1980_2016_DM.nc

# Regrid, set missing value to zero, and adjust attributes
cdo -P 6 -f nc4   setmisstoc,0 -setattribute,climatic_mass_balance@units="kg m-2 day-1" -chname,gld,climatic_mass_balance -chname,tas,ice_surface_temp -chname,pr,precipitation -remapycon,${grid_file} -setgrid,../../grids/rotated_grid.txt DMI-HIRHAM5_GL2_ERAI_1980_2016_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_EPSG3413_${grid}M_DM.nc

