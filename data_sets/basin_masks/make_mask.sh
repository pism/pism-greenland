#!/bin/bash

UGID=225
ores=$1

x_min=-282650.0
x_max=293350.0
y_min=-2417600.0
y_max=-2021600.0

gdal_rasterize  -where "UGID='${UGID}'" -l Greenland_Basins_PS_v1.4.2_1980_epsg3413 -tr $ores $ores -a_nodata -9999 -burn 0 -te $x_min $y_min $x_max $y_max -ot Byte ~/base/gris-analysis/basins/Greenland_Basins_PS_v1.4.2_1980_epsg3413.shp ugid_225_Jakobshavn_Isbrae_mask_epsg3413_g${ores}m.nc

