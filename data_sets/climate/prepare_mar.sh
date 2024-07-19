#!/bin/bash

mar_dir=mar_v3.14_era

mkdir -p $mar_dir

start_year=1975
end_year=2023
for year in $(seq $start_year $end_year); do
    wget -nc --no-check-certificate -P $mar_dir http://ftp.climato.be/fettweis/MARv3.14/Greenland/ERA5-1km-monthly/MARv3.14-monthly-ERA5-${year}.nc
done

exit 
cdo -f nc4 -z zip_2 chname,SMBcorr,climatic_mass_balance,STcorr,ice_surface_temp,T2Mcorr,air_temp,RUcorr,water_input_rate -setattribute,SMBcorr@units="kg m-2 month-1",RUcorr@units="kg m-2 month-1" -selvar,SMBcorr,T2Mcorr,STcorr,RUcorr -mergetime $mar_dir/MARv3.14-monthly-ERA5-*.nc MARv3.14-monthly-ERA5-${start_year}_${end_year}.nc

adjust_timeline.py -p monthly -a ${start_year}-01-01 -d ${end_year}-01-01 MARv3.14-monthly-ERA5-${start_year}_${end_year}.nc

exit

cdo -f nc4 -z zip_2 chname,SMBcorr,climatic_mass_balance,STcorr,ice_surface_temp,TTcorr,air_temp,RUcorr,water_input_rate -setattribute,SMBcorr@units="kg m-2 month-1",RUcorr@units="kg m-2 month-1" -selvar,SMBcorr,TTcorr,STcorr,RUcorr -mergetime 20CRv2c_1900-2014_20km/MARv3.5.2-20km-monthly-20CRv2c-*.nc MARv3.5.2-20km-monthly-20CRv2c_1900_2014.nc

adjust_timeline.py -p monthly -a 1900-01-01 -d 1900-01-01 MARv3.5.2-20km-monthly-20CRv2c_1900_2014.nc

cdo -O -f nc4 -z zip_2 setmisstoc,0 -remapbil,../../grids/gris_ext_g4500m.txt -setgrid,../../grids/mar_grid.txt  MARv3.5.2-20km-monthly-20CRv2c_1900_2014.nc  MARv3.5.2-20km-monthly-20CRv2c_4500m_1900_2014.nc
