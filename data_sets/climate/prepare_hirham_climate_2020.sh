#!/bin/bash

hh5_dir=hirham5_zip
hh5_nc_dir=hirham5_nc
mkdir -p $hh5_dir
mkdir -p $hh5_nc_dir

start_year=1980
end_year=2021
for year in $(seq $start_year $end_year); do
    wget -nc --no-check-certificate -P $hh5_dir http://ensemblesrt3.dmi.dk/data/prudence/temp/nichan/Daily2D_GrIS/${year}.zip
    unzip -u ${hh5_dir}/${year}.zip
    cdo -O -f nc4 -z zip_2 -r settbounds,day -setattribute,gld@units="kg m-2 day-1",rfrz@units="kg m-2 day-1",snmel@units="kg m-2 day-1",rogl@units="kg m-2 day-1",snfall@units="kg m-2 day-1",rainfall@units="kg m-2 day-1" -selvar,gld,tas,rogl,snfall,rainfall,rfrz,sn,snmel -setmissval,1e19 -setgrid,../../grids/rotated_grid.txt -mergetime ${year}/Daily2D_*.nc ${hh5_nc_dir}/${year}.nc
    rm -rf ${year}/Daily2D_*.nc
done

cdo -O -f nc4 -z zip_2 mergetime ${hh5_nc_dir}/*.nc DMI-HIRHAM5_${start_year}_${end_year}.nc

cdo -O -P 6 -f nc4  selvar,climatic_mass_balance,air_temp,ice_surface_temp,precipitation,water_input_rate -setattribute,precipitation@units="kg m-2 day-1" -aexpr,"precipitation=snfall+rainfall;air_temp=ice_surface_temp" -chname,rogl,water_input_rate -chname,gld,climatic_mass_balance -chname,tas,ice_surface_temp -remapycon,../../grids/gris_ext_g4500m.txt -setmisstodis -setgrid,../../grids/rotated_grid.txt DMI-HIRHAM5_${start_year}_${end_year}.nc DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM.nc

cdo -f nc4 -z zip_2 -O selyear,1980/1984 DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM.nc DMI-HIRHAM5_ERA_1975_1979_EPSG3413_4500M_DM.nc
~/pism/sources/util/adjust_timeline.py -p daily -a 1975-01-01 -d ${start_year}-01-01 DMI-HIRHAM5_ERA_1975_1979_EPSG3413_4500M_DM.nc

cdo -f nc4 -z zip_2 -O mergetime -selyear,1975/1979 DMI-HIRHAM5_ERA_1975_1979_EPSG3413_4500M_DM.nc  DMI-HIRHAM5_ERA_${start_year}_${end_year}_EPSG3413_4500M_DM.nc  DMI-HIRHAM5_ERA_1975_${end_year}_EPSG3413_4500M_DM.nc

cdo -O -f nc4 -z zip_2 monmean DMI-HIRHAM5_ERA_1975_${end_year}_EPSG3413_4500M_DM.nc DMI-HIRHAM5_ERA_1975_${end_year}_EPSG3413_4500M_MM.nc
cdo -O -f nc4 -z zip_2 ymonmean DMI-HIRHAM5_ERA_1975_${end_year}_EPSG3413_4500M_DM.nc DMI-HIRHAM5_ERA_1975_${end_year}_EPSG3413_4500M_YMM.nc
cdo -O -f nc4 -z zip_2 timmean DMI-HIRHAM5_ERA_1975_${end_year}_EPSG3413_4500M_DM.nc DMI-HIRHAM5_ERA_1975_${end_year}_EPSG3413_4500M_TM.nc
