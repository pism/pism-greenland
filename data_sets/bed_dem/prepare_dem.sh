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


run_with_mpi () {

    if [ -z "$SLURM_JOBID" ];
    then
        echo "running without SLURM $*"
        mpiexec -n $*
    else
        echo "running under SLURM $*"
        mpirun -machinefile ./nodes_$SLURM_JOBID -np $*
    fi
}


# First we download the Bamber 2001 SeaRISE data set

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc
echo "fetching master file ... "
wget -nc ${DATAURL}${DATANAME}   # -nc is "no clobber"
echo "  ... done."
echo
PISMVERSION=pism_$DATANAME
echo -n "creating bootstrapable $PISMVERSION from $DATANAME ... "
# copy the vars we want, and preserve history and global attrs
ncks -O -v mapping,lat,lon,bheatflx,topg,thk,presprcp,smb,airtemp2m $DATANAME $PISMVERSION
# convert from m yr-1 to kg m-2 yr-1
ncap2 -O -s "precipitation=presprcp*1000.0" $PISMVERSION $PISMVERSION
ncatted -O -a units,precipitation,o,c,"kg m-2 yr-1" $PISMVERSION
ncatted -O -a long_name,precipitation,c,c,"mean annual precipitation rate" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION
# use pism-recognized name for 2m air temp
ncrename -O -v airtemp2m,ice_surface_temp  $PISMVERSION
ncatted -O -a units,ice_surface_temp,c,c,"Celsius" $PISMVERSION
# use pism-recognized name and standard_name for surface mass balance, after
# converting from liquid water equivalent thickness per year to [kg m-2 year-1]
ncap2 -t $NN -O -s "climatic_mass_balance=1000.0*smb" $PISMVERSION $PISMVERSION
ncatted -O -a standard_name,climatic_mass_balance,m,c,"land_ice_surface_specific_mass_balance" $PISMVERSION
ncatted -O -a units,climatic_mass_balance,m,c,"kg m-2 year-1" $PISMVERSION
# This is a *choice* of the model of surface mass balance in thk==0 areas.
ncap2 -O -s "where(thk <= 0.0){climatic_mass_balance=-10000.0;}" $PISMVERSION $PISMVERSION
# de-clutter by only keeping vars we want
ncks -O -v mapping,lat,lon,bheatflx,topg,thk,precipitation,ice_surface_temp,climatic_mass_balance \
  $PISMVERSION $PISMVERSION
# straighten dimension names
ncrename -O -d x1,x -d y1,y -v x1,x -v y1,y $PISMVERSION $PISMVERSION
nc2cdo.py $PISMVERSION
echo "done."
echo


ibcaofile=IBCAO_V3_500m_RR
wget -nc http://www.ngdc.noaa.gov/mgg/bathymetry/arctic/grids/version3_0/${ibcaofile}_tif.zip
#unzip -o ${ibcaofile}_tif.zip

ismip6_grid() {

# buffer = grid spacing / 2
buffer_x=500
buffer_y=500
xmin=$((-720000 - $buffer_x))
ymin=$((-3450000 - $buffer_y))
xmax=$((960000 + $buffer_x))
ymax=$((-570000 + $buffer_y))


for GRID in 1000; do
    outfile_prefix=pism_Greenland_ismip6_${GRID}m_v${ver}
    outfile=${outfile_prefix}.nc
    outfile_ctrl=${outfile_prefix}_ctrl.nc
    outfile_nb=${outfile_prefix}_wc.nc
    outfile_rm=${outfile_prefix}_rm.nc
    outfile_sm_prefix=pism_Greenland_${GRID}m_v${ver}
    outfile_sm_ctrl=${outfile_sm_prefix}_ctrl.nc
    outfile_sm_nb=${outfile_sm_prefix}_wc.nc
    outfile_sm_rm=${outfile_sm_prefix}_rm.nc
    
    for var in "bed"; do
        rm -f g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc
        gdalwarp $CUT -overwrite  -r average -s_srs EPSG:3413 -t_srs EPSG:3413 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID -of GTiff NETCDF:$infile:$var g${GRID}m_${var}_v${ver}.tif
        gdal_translate -a_srs epsg:3413 -co "FORMAT=NC4" -of netCDF g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc 
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
        gdal_translate -a_srs epsg:3413 -co "FORMAT=NC4" -of netCDF g${GRID}m_${var}_v${ver}.tif g${GRID}m_${var}_v${ver}.nc 
    done
    
    ncks -O g${GRID}m_bed_v${ver}.nc $outfile
    ncatted -a _FillValue,bed,d,, $outfile
    for var in "surface" "thickness" "mask"; do
        ncks -A g${GRID}m_${var}_v${ver}.nc $outfile
    done
        
    # This is not needed, but it can be used by PISM to calculate correct cell volumes, and for remapping scripts"
    ncatted -a proj,global,o,c,"epsg:3413" $outfile

    # Add IBCAO bathymetry for the outer part of the domain
    gdalwarp $CUT -overwrite -r average -t_srs EPSG:3413 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID -of GTiff ${ibcaofile}_tif/${ibcaofile}.tif ${ibcaofile}_epsg3413_g${GRID}m.tif
    gdal_translate -co "FORMAT=NC4" -of netCDF  ${ibcaofile}_epsg3413_g${GRID}m.tif  ${ibcaofile}_epsg3413_g${GRID}m.nc
    ncks -4 -A -v Band1 ${ibcaofile}_epsg3413_g${GRID}m.nc $outfile
    ncap2 -O -s "where(bed==-9999) {bed=Band1;}; where(Band1<=-9990) {bed=-9999;};" $outfile $outfile
    ncks -O -v Band1,topg -x $outfile $outfile

    # Remap SeaRISE fields
    ncks -4 -O g${GRID}m_${var}_v${ver}.nc griddes_${GRID}m.nc
    nc2cdo.py --srs "epsg:3413" griddes_${GRID}m.nc
    if [[ $N == 1 ]] ; then
        cdo -f nc4 remapbil,griddes_${GRID}m.nc ${PISMVERSION} v${ver}_tmp_${GRID}m_searise.nc
    else
        cdo -P $N -f nc4 remapbil,griddes_${GRID}m.nc ${PISMVERSION} v${ver}_tmp_${GRID}m_searise.nc
    fi

    # Extrapolate climate and other BC fields onto PISM domain
    run_with_mpi $NN fill_missing_petsc.py -v precipitation,ice_surface_temp,bheatflx,climatic_mass_balance v${ver}_tmp_${GRID}m_searise.nc v${ver}_tmp2_${GRID}m.nc
    ncap2 -O -s "polar_stereographic=char(polar_stereographic);"  v${ver}_tmp2_${GRID}m.nc  v${ver}_tmp2_${GRID}m.nc
    ncks -4 -A -C -v precipitation,ice_surface_temp,bheatflx,climatic_mass_balance v${ver}_tmp2_${GRID}m.nc $outfile
    ncatted -a _FillValue,,d,, -a missing_value,,d,, $outfile
    
    # remove regridding artifacts, give precedence to mask: we set thickness and
    # surface to 0 where mask has ocean
    # Make FTT mask: 1 where there is no floating (3) or grounded (2) ice

    # Instead of Ocean Kill we now use "land area fraction"

    ncap2 -O -s "where(mask==2) thickness=surface-bed; where(thickness<0) thickness=0; ftt_mask[\$y,\$x]=0b; where(mask==0) {thickness=0.; surface=0.;}; where(mask!=2) ftt_mask=1; where(mask!=3) ftt_mask=1;" $outfile $outfile
    ncap2 -O -s 'land_ice_area_fraction_retreat = thickness; where(thickness > 0 || thickness + bed >= (1 - 910.0/1028.0) * thickness + 0) land_ice_area_fraction_retreat = 1;land_ice_area_fraction_retreat@units="1";land_ice_area_fraction_retreat@long_name="maximum ice extent mask";land_ice_area_fraction_retreat@standard_name="";' $outfile $outfile

    ncks -h -O $outfile $outfile_ctrl
    ncks -h -O $outfile $outfile_nb

    # Here we cut out the topography of the Canadian Archipelago
    var=thickness
    gdalwarp -overwrite -dstnodata 0 -cutline  ../shape_files/gris-domain-ismip6.shp NETCDF:$outfile_nb:$var g${GRID}m_nb_${var}_v${ver}.tif
    gdal_translate -of netCDF -co "FORMAT=NC4" g${GRID}m_nb_${var}_v${ver}.tif g${GRID}m_nb_${var}_v${ver}.nc
    ncks -A -v $var g${GRID}m_nb_${var}_v${ver}.nc $outfile_nb
    var=bed
    gdalwarp -overwrite -dstnodata -9999 -cutline  ../shape_files/gris-domain-ismip6.shp NETCDF:$outfile_nb:$var g${GRID}m_nb_${var}_v${ver}.tif
    gdal_translate -of netCDF -co "FORMAT=NC4" g${GRID}m_nb_${var}_v${ver}.tif g${GRID}m_nb_${var}_v${ver}.nc
    ncks -A -C -v $var g${GRID}m_nb_${var}_v${ver}.nc $outfile_nb
    ncatted -a _FillValue,bed,d,, -a _FillValue,thickness,d,, $outfile_nb
    ncap2 -O -s "where(bed==-9999) {mask=0; surface=0; thickness=0;};"  $outfile_nb  $outfile_nb
    nccopy $outfile_nc $outfile_rm

done

}

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
    
    ncks -O g${GRID}m_bed_v${ver}.nc $outfile
    ncatted -a _FillValue,bed,d,, $outfile
    for var in "surface" "thickness" "mask"; do
        ncks -A g${GRID}m_${var}_v${ver}.nc $outfile
    done
        
    # This is not needed, but it can be used by PISM to calculate correct cell volumes, and for remapping scripts"
    ncatted -a proj,global,o,c,"epsg:3413" $outfile

    # Add IBCAO bathymetry for the outer part of the domain
    gdalwarp $CUT -overwrite -r average -t_srs EPSG:3413 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID -of GTiff ${ibcaofile}_tif/${ibcaofile}.tif ${ibcaofile}_epsg3413_g${GRID}m.tif
    gdal_translate -co "FORMAT=NC4" -of netCDF  ${ibcaofile}_epsg3413_g${GRID}m.tif  ${ibcaofile}_epsg3413_g${GRID}m.nc
    ncks -A -v Band1 ${ibcaofile}_epsg3413_g${GRID}m.nc $outfile
    ncap2 -O -s "where(bed==-9999) {bed=Band1;}; where(Band1<=-9990) {bed=-9999;};" $outfile $outfile
    ncks -O -v Band1,topg -x $outfile $outfile

    # Remap SeaRISE fields
    ncks -4 -O g${GRID}m_mask_v${ver}.nc griddes_${GRID}m.nc
    nc2cdo.py --srs "epsg:3413" griddes_${GRID}m.nc
    if [[ $N == 1 ]] ; then
        cdo -f nc4 remapbil,griddes_${GRID}m.nc ${PISMVERSION} v${ver}_tmp_${GRID}m_searise.nc
    else
        cdo -P $N -f nc4 remapbil,griddes_${GRID}m.nc ${PISMVERSION} v${ver}_tmp_${GRID}m_searise.nc
    fi

    # Extrapolate climate and other BC fields onto PISM domain
    run_with_mpi $NN fill_missing_petsc.py -v precipitation,ice_surface_temp,bheatflx,climatic_mass_balance v${ver}_tmp_${GRID}m_searise.nc v${ver}_tmp2_${GRID}m.nc
    ncap2 -O -s "polar_stereographic=char(polar_stereographic);"  v${ver}_tmp2_${GRID}m.nc  v${ver}_tmp2_${GRID}m.nc
    ncks -4 -A -v precipitation,ice_surface_temp,bheatflx,climatic_mass_balance v${ver}_tmp2_${GRID}m.nc $outfile
    ncatted -a _FillValue,,d,, -a missing_value,,d,, $outfile
    ncks -A -v time Greenland_5km_v1.1.nc  $outfile
    
    # remove regridding artifacts, give precedence to mask: we set thickness and
    # surface to 0 where mask has ocean
    # Make FTT mask: 1 where there is no floating (3) or grounded (2) ice

    # Instead of Ocean Kill we now use "land area fraction"
    ncap2 -O -s "where(mask==2) thickness=surface-bed; where(thickness<0) thickness=0; ftt_mask[\$y,\$x]=0b; where(mask==0) {thickness=0.; surface=0.;}; where(mask!=2) ftt_mask=1; where(mask!=3) ftt_mask=1;" $outfile $outfile
    ncap2 -O -s 'land_ice_area_fraction_retreat = thickness; where(thickness > 0 || thickness + bed >= (1 - 910.0/1028.0) * thickness + 0) land_ice_area_fraction_retreat = 1;land_ice_area_fraction_retreat@units="1";land_ice_area_fraction_retreat@long_name="maximum ice extent mask";land_ice_area_fraction_retreat@standard_name="";' $outfile $outfile


    ncks -h -O $outfile $outfile_ctrl
    ncks -h -O $outfile $outfile_nb

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
    ncap2 -O -s "where(bed==-9999) {mask=0; surface=0; thickness=0;};"  $outfile_nb  $outfile_nb

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

done

}

# ismip6_grid
default_grid
