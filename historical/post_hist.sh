#!/bin/bash

grid=900
model=HH5-MIROC85-YM
for d in jib nw; do
cdo -f nc4 -z zip_3 yearmean 2019_04_hindcast_tmp/ex_${d}_g${grid}m_v3a_id_${model}_2008-1-1_2017-1-1.nc  2019_04_hindcast/spatial/ex_${d}_g${grid}m_v3a_id_${model}_2008-1-1_2017-1-1_YM.nc
extract_interface.py -t ice_ocean -o 2019_04_hindcast/ice_ocean/io_${d}_g${grid}m_v3a_id_${model}_2008-1-1_2017-1-1_YM.nc 2019_04_hindcast/spatial/ex_${d}_g${grid}m_v3a_id_${model}_2008-1-1_2017-1-1_YM.nc
done
