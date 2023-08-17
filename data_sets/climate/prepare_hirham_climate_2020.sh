#!/bin/bash

hh5_dir=hirham5_zip
hh5_nc_dir=hirham5_nc
mkdir -p $hh5_dir
mkdir -p $hh5_nc_dir

start_year=1980
end_year=2020
for year in $(seq $start_year $end_year); do
    wget -nc -P $hh5_dir http://ensemblesrt3.dmi.dk/data/prudence/temp/nichan/Daily2D_GrIS/${year}.zip
    # unzip -u ${hh5_dir}/${year}.zip
    cdo -O -f nc4 -z zip_2 -r settbounds,day -setattribute,gld@units="kg m-2 day-1",rfrz@units="kg m-2 day-1",snmel@units="kg m-2 day-1",rogl@units="kg m-2 day-1",snfall@units="kg m-2 day-1",rainfall@units="kg m-2 day-1" -selvar,gld,tas,rogl,snfall,rainfall,rfrz,sn,snmel -setmissval,1e19 -setgrid,../../grids/rotated_grid.txt -mergetime ${year}/Daily2D_*.nc ${hh5_nc_dir}/${year}.nc
done
cdo -O -f nc4 -z zip_2 mergetime ${hh5_nc_dir}/*.nc DMI-HIRHAM5_${start_year}_${end_year}.nc
cdo -f nc4 -z zip_2 ydaymean DMI-HIRHAM5_${start_year}_${end_year}.nc  DMI-HIRHAM5_${start_year}_${end_year}_YDM.nc


cdo -O -P 6 -f nc4  selvar,climatic_mass_balance,air_temp,ice_surface_temp,precipitation,water_input_rate -setattribute,precipitation@units="kg m-2 day-1" -aexpr,"precipitation=snfall+rainfall;air_temp=ice_surface_temp" -chname,rogl,water_input_rate -chname,gld,climatic_mass_balance -chname,tas,ice_surface_temp -remapbil,../../grids/gris_ext_g4500m.txt -setgrid,../../grids/rotated_grid.txt  DMI-HIRHAM5_${start_year}_${end_year}.nc DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM_all.nc

cdo -O  -f nc4  -z zip_2 setmisstoc,0 -selvar,climatic_mass_balance,precipitation,water_input_rate  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM_all.nc  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM_rest.nc
cdo -O  -f nc4  -z zip_2 setmisstoc,270 -selvar,ice_surface_temp,air_temp  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM_all.nc  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM_temp.nc

cdo -O  -f nc4 -z zip_2 merge  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM_rest.nc  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM_temp.nc  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM.nc
cdo -O -f nc4 -z zip_2 monmean DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM.nc DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_MM.nc
cdo -O -f nc4 -z zip_2 ymonmean DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM.nc DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_YMM.nc
cdo -O -f nc4 -z zip_2 timmean DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM.nc DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_TM.nc
