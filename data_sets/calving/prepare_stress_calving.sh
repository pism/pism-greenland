#!/bin/bash

# csvsql --db sqlite:///myjoindb.db --insert myjoincsv.csv
# ogr2ogr -append -f "SQLite" myjoindb.db myjoinshp.shp
# ogr2ogr -f "ESRI Shapefile" -sql "SELECT csv.*, shp.* FROM myjoinshp shp INNER JOIN myjoincsv csv ON csv.joinfield = shp.joinfield" joined_output.shp myjoindb.db

basin_file=$1
vc_file=$2

csvsql --db sqlite:///vonmises_calving_threshold.db --insert "$vc_file"
ogr2ogr -append -f "SQLite" vonmises_calving_threshold.db "$basin_file"

ogr2ogr -f "ESRI Shapefile" -sql "SELECT csv.*, shp.* FROM "$basin_file" shp INNER JOIN ugid_vc csv ON csv.UGID = shp.UGID" vonmises_calving_threshold.shp vonmises_calving_threshold.db
