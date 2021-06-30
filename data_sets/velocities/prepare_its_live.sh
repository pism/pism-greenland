#!/bin/bash

for year in {1985..2018}; do
    echo "Preparing year $year"
    cdo -L -f nc4 settbounds,year -setdate,${year}-1-1 -setreftime,${year}-1-1 itslive/GRE_G0240_${year}.nc GRE_G0240_T_${year}.nc 
done

file=GRE_G0240_1985_2018.nc
rm -f $file

cdo -f nc4 -z zip_2 chname,v,velsurf_mag -mergetime GRE_G0240_*.nc $file

extract_profiles.py -v velsurf_mag --srs epsg:3413 /Volumes/zachariae/crios2pism/data_sets/shape_files/joughin-gps-points.shp $file profiles/gps_stations_$file
