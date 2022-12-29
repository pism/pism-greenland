#!/bin/bash

# (c) 2015-2019 Andy Aschwanden

# The script prepares the Greenland BedMachine for use with
# the Parallel Ice Sheet Model (PISM)
# The user needs to download the BedMachine file from NSIDC
# https://nsidc.org/data/idbmg4
# e.g. run with
# sh preprocess_mc_bed.sh 4

# Update Source for DEM for RAGIS

set -x -e

export HDF5_USE_FILE_LOCKING=FALSE

# run ./preprocess.sh 1 if you haven't CDO compiled with OpenMP
NN=4  # default number of processors
if [ $# -gt 0 ] ; then
  NN="$1"
fi
N=$NN


infile=BedMachineGreenland-v5.nc
if [ -n "$2" ]; then
    infile=$2
fi

ver=2022
if [ -n "$3" ]; then
    ver=$3
fi


gebco=GEBCO_2022_sub_ice_topo.nc



default_grid() {

# Create a buffer that is a multiple of the grid resolution
# and works for grid resolutions up to 36km.
buffer_x=148650
buffer_y=130000
xmin=$((-638000 - $buffer_x - 468000))
ymin=$((-3349600 - $buffer_y))
xmax=$((864700 + $buffer_x))
ymax=$((-657600 + $buffer_y))

for GRID in 18000 9000 6000 4500 3600 3000 2400 1800 1500 1200 900 600 450 300 150; do
    outfile_prefix=pism_Greenland_ext_${GRID}m_mcb_jpl_v${ver}
    outfile=${outfile_prefix}.nc
    outfile_ctrl=${outfile_prefix}_ctrl.nc
    outfile_nb=${outfile_prefix}_wc.nc
    outfile_sm_prefix=pism_Greenland_${GRID}m_mcb_jpl_v${ver}
    outfile_sm_ctrl=${outfile_sm_prefix}_ctrl.nc
    outfile_sm_nb=${outfile_sm_prefix}_wc.nc
    
    for var in "bed"; do
        rm -f g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc
        gdalwarp $CUT -overwrite  -r average -s_srs EPSG:3413 -t_srs EPSG:3413 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID -of GTiff NETCDF:$infile:$var g${GRID}m_${var}_v${ver}.tif
        gdal_translate -co "FORMAT=NC4" -of netCDF g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc 
        ncatted -a nx,global,d,, -a ny,global,d,, -a xmin,global,d,, -a ymax,global,d,, -a spacing,global,d,, g${GRID}m_${var}_v${ver}.nc
        
    done
    for var in "surface" "thickness"; do
        rm -f g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc
        gdalwarp -overwrite -r average -te $xmin $ymin $xmax $ymax -tr $GRID $GRID -of GTiff NETCDF:$infile:$var g${GRID}m_${var}_v${ver}.tif
        gdal_translate -a_srs epsg:3413 -co "FORMAT=NC4" -of netCDF g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc
        ncatted -a _FillValue,$var,d,, g${GRID}m_${var}_v${ver}.nc
        ncap2 -O -s "where(${var}<=0) ${var}=0.;" g${GRID}m_${var}_v${ver}.nc g${GRID}m_${var}_v${ver}.nc
    done
    for var in "mask"; do
        rm -f g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc
        gdalwarp -overwrite -r near -te $xmin $ymin $xmax $ymax -tr $GRID $GRID -of GTiff NETCDF:$infile:$var g${GRID}m_${var}_v${ver}.tif
        gdal_translate -co "FORMAT=NC4" -of netCDF g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc 
    done

    for var in "elevation"; do
        rm -f g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc
        gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:3413 -r average -te $xmin $ymin $xmax $ymax -tr $GRID $GRID -of GTiff NETCDF:$gebco:$var g${GRID}m_${var}_v${ver}.tif
        gdal_translate -co "FORMAT=NC4" -of netCDF g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc 
    done

    ncks -O g${GRID}m_bed_v${ver}.nc $outfile
    ncatted -a _FillValue,bed,d,, $outfile
    for var in "surface" "thickness" "mask" "elevation"; do
        ncks -A g${GRID}m_${var}_v${ver}.nc $outfile
    done
    ncap2 -O -s "where(bed==-9999) bed=elevation;" $outfile $outfile 
    # This is not needed, but it can be used by PISM to calculate correct cell volumes, and for remapping scripts"
    ncatted -a proj,global,o,c,"epsg:3413" $outfile

    # Instead of Ocean Kill we now use "land area fraction"
    ncap2 -O -s "where(mask==2) thickness=surface-bed; where(thickness<0) thickness=0; ftt_mask[\$y,\$x]=0b; where(mask==0) {thickness=0.; surface=0.;}; where(mask!=2) ftt_mask=1; where(mask!=3) ftt_mask=1;" $outfile $outfile
    ncap2 -O -s 'land_ice_area_fraction_retreat = thickness; where(thickness > 0 || thickness + bed >= (1 - 910.0/1028.0) * thickness + 0) land_ice_area_fraction_retreat = 1;land_ice_area_fraction_retreat@units="1";land_ice_area_fraction_retreat@long_name="maximum ice extent mask";land_ice_area_fraction_retreat@standard_name="";' $outfile $outfile

    ncks -h -O -4 -L 2 -v elevation -x $outfile $outfile
    ncks -h -O  $outfile $outfile_ctrl
    ncks -h -O  $outfile $outfile_nb

    ncks -O -4 -L 2 -v elevation -x $outfile $outfile
    
    # Here we cut out the topography of the Canadian Archipelago
    var=thickness
    gdalwarp -overwrite -dstnodata 0 -cutline  ../shape_files/gris-domain-paleo.shp NETCDF:$outfile_nb:$var g${GRID}m_nb_${var}_v${ver}.tif
    gdal_translate -of netCDF -co "FORMAT=NC4" g${GRID}m_nb_${var}_v${ver}.tif g${GRID}m_nb_${var}_v${ver}.nc
    ncks -A -v $var g${GRID}m_nb_${var}_v${ver}.nc $outfile_nb
    var=bed
    gdalwarp -overwrite -dstnodata -9999 -cutline  ../shape_files/gris-domain-paleo.shp NETCDF:$outfile_nb:$var g${GRID}m_nb_${var}_v${ver}.tif
    gdal_translate -of netCDF -co "FORMAT=NC4" g${GRID}m_nb_${var}_v${ver}.tif g${GRID}m_nb_${var}_v${ver}.nc
    ncks -A -C -v $var g${GRID}m_nb_${var}_v${ver}.nc $outfile_nb
    ncatted -a _FillValue,bed,d,, -a _FillValue,thickness,d,, $outfile_nb
    ncap2 -O -s "where(bed==-9999) {mask=0; surface=0; thickness=0; bed=-300;};"  $outfile_nb  $outfile_nb

    
    # Cut out smaller domain used for projections
    e0=-638000
    n0=-3349600
    e1=864700
    n1=-657600

    buffer_e=40650
    buffer_n=22000
    e0=$(($e0 - $buffer_e))
    n0=$(($n0 - $buffer_n))
    e1=$(($e1 + $buffer_e))
    n1=$(($n1 + $buffer_n))

    # Shift to cell centers
    e0=$(($e0 + $GRID / 2 ))
    n0=$(($n0 + $GRID / 2))
    e1=$(($e1 - $GRID / 2))
    n1=$(($n1 - $GRID / 2))

    ncks -O -d x,$e0.,$e1. -d y,$n0.,$n1.  $outfile_ctrl  $outfile_sm_ctrl
    ncks -O -d x,$e0.,$e1. -d y,$n0.,$n1.  $outfile_nb  $outfile_sm_nb

    # rm -f  g${GRID}m_*_v${ver}.nc
   
done

}

# ismip6_grid
default_grid
