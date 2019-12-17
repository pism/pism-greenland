#!/bin/bash

cdo setattribute,water_input_rate@units="Gt yr-1" -divc,1e12 -fldsum -mulc,1e6 -seltimestep,1/365 -selvar,water_input_rate synth_jib_runoff_g1000m.nc fldsum_synth_jib_runoff_g1000m.nc 
for m in fldmax fldmean fldsum; do
cdo -f nc4 -z zip_3 ${m} -selyear,10 2019_12_routing_tmp/ex_synth_jib_g1000m_id_ROUTING_0_10.nc 2019_12_routing/spatial/${m}_ex_synth_jib_g1000m_id_ROUTING_0_10.nc
cdo -f nc4 -z zip_3 ${m} -selyear,10 2019_12_steady_tmp/ex_synth_jib_g1000m_id_STEADY_0_10.nc 2019_12_steady/spatial/${m}_ex_synth_jib_g1000m_id_STEADY_0_10.nc
done
