#!/bin/bash

x_min=-280000.0
x_max=320000.0
y_min=-2410000.0
y_max=-2020000.0

ncks -d x,$x_min,$x_max -d y,$y_min,$y_max pism_Greenland_4500m_mcb_jpl_v3a_ctrl.nc pism_JIB_4500m_mcb_jpl_v3a_ctrl.nc

ncks -d x,$x_min,$x_max -d y,$y_min,$y_max DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_YMM_EPSG3413_3600m_0.nc JIB_DMI-HIRHAM5_GL2_ERAI_1980_2014_MRROS_YMM_EPSG3413_3600m_0.nc
