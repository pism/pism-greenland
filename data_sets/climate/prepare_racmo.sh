#!/bin/bash

# Download files from https://surfdrive.surf.nl/files/index.php/s/No8LoNA18eS1v72


racmo_dir=Monthly

cdo -O -f nc4 -z zip_2 selyear,1975/2023 -selvar,t2m,precip,smb,runoff -merge $racmo_dir/precip_monthlyS_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202303.nc $racmo_dir/smb_monthlyS_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202303.nc $racmo_dir/runoff_monthlyS_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202303.nc $racmo_dir/t2m_monthlyA_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202303.nc monthly_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202303.nc


cdo -O -P 6 -f nc4 remapycon,../../grids/gris_ext_g4500m.txt -setmisstodis -setgrid,../../grids/rotated_grid_racmo.txt monthly_FGRN055_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202303.nc RACMO2.3p2_ERA5_1975_2023.nc
