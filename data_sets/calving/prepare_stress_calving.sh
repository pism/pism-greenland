#!/bin/bash

basin_file=$1
csvsql --db sqlite:///vonmises_calving_threshold.db --insert ugid_vc.csv
ogr2ogr -append -f "SQLite" vonmises_calving_threshold.db $basin_file

ogr2ogr -f "ESRI Shapefile" -sql "SELECT csv.*, shp.* FROM myjoinshp shp INNER JOIN ugid_vc csv ON csv.joinfield = shp.joinfield" vonmises_calving_threshold.shp vonmises_calving_threshold.db
