#!/bin/bash

for var in PR TAS; do
    for file in DMI-HIRHAM5_GL2_ERAI_1980_2014_${var}_DM.nc DMI-HIRHAM5_GL2_ERAI_2015_2016_${var}_DM.nc; do 
        echo "Downloading $file"
        wget -nc http://prudence.dmi.dk/data/temp/RUM/HIRHAM/GREENLAND/ERAI/DAILY/$file.gz
        echo "Unzipping ${file}.gz"
        gunzip ${file}.gz
        echo "Removing height dimension"
        ncwa -O -a height $file $file
    done
    cdo --reduce_dim sellevel,0 -mergetime DMI-HIRHAM5_GL2_ERAI_1980_2014_${var}_DM.nc DMI-HIRHAM5_GL2_ERAI_2015_2016_${var}_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2014_${var}_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_${var}_DM.nc
done

wget -nc http://prudence.dmi.dk/data/temp/RUM/HIRHAM/GREENLAND/ERAI/DAILY/DMI-HIRHAM5_GL2_ERAI_1980_2016_gld_DM.nc.gz
gunzip DMI-HIRHAM5_GL2_ERAI_1980_2016_gld_DM.nc.gz

cdo -O -f nc4 merge  DMI-HIRHAM5_GL2_ERAI_1980_2016_gld_DM.nc  DMI-HIRHAM5_GL2_ERAI_1980_2016_TAS_DM.nc  DMI-HIRHAM5_GL2_ERAI_1980_2016_PR_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_DM.nc

create_greenland_epsg3413_grid.py -g 4500 g4500m.nc
nc2cdo.py g4500m.nc

adjust_timeline.py -p daily -a 1980-1-1 -d 1980-1-1 DMI-HIRHAM5_GL2_ERAI_1980_2016_DM.nc 
cdo -f nc4 -L  setmisstoc,0 -setattribute,climatic_mass_balance@units="kg m-2 day-1" -chname,gld,climatic_mass_balance -chname,tas,ice_surface_temp -chname,pr,precipitation -remapycon,g4500m.nc -setgrid,b../../util/rotated_grid.txt DMI-HIRHAM5_GL2_ERAI_1980_2016_DM.nc DMI-HIRHAM5_GL2_ERAI_1980_2016_EPSG3413_4500M_DM.nc
