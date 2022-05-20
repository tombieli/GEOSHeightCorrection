#!/bin/bash

error_handle()
{
        echo [`date`] Unspecified error in line $1 >&2
        echo Exiting >&2
        if [ -e $TEMP_DIR ];
        then
                $RM -rf $TEMP_DIR
        fi
        exit 2
}



trap 'error_handle ${LINENO}' ERR

getULLR()
{
        local IN_LINES=`$GDAL_INFO $1 | grep "Upper Left  \|Lower Right " | sed -E 's#^[^0-9]*([1-9][0-9]*[.]?[0-9]+), ([1-9][0-9]*[.]?[0-9]+).*$#\1 \2#'`
        for L in $IN_LINES; do
                echo "$L "
        done
}

getPixelSize()
{
	$GDAL_INFO $1 | grep "Pixel Size " | sed -E 's#^[^0-9]*([1-9][0-9]*[.]?[0-9]+),[-]?([1-9][0-9]*[.]?[0-9]+).*$#\1 \2#'
}

getSize()
{
	$GDAL_INFO $1 | grep "Size is " | sed -E 's#^[^0-9]*([1-9][0-9]*), ([1-9][0-9]*).*$#\1 \2#'
}


getBandsCount()
{
	$GDAL_INFO $1 | grep "Band" | cut -f 2 -d " " | tail -1
}

getNoData()
{
	$GDAL_INFO $1 | grep NoData | head -$2 | tail -1 | sed -E 's#.*NoData Value=(.*)$#\1#'
}

getNoDataTest()
{
	local V=$1
	local NO_DATA_V=$2
	if [ "$NO_DATA_V" == nan ];
	then
		echo $V=$V
	else
		echo $V!=${NO_DATA_V}
	fi
}

getProjSRS()
{
	$GDAL_INFO -proj4 $1 | grep -A 1 "PROJ" | tail -1 | tr -d "'"
}

getXYZAllBandsOption()
{
	GV=$[${GDAL_VERSION[0]} * 10 + ${GDAL_VERSION[1]}]
	if [ "$GV" -ge 33 ]; then
		echo "-allbands"
	else
		seq $1 | sed "s#^#-band #"
	fi
}

errcheck()
{
	echo $@ >&2
        if ! $@; then
                echo "\"$@\"" failed >&2
                exit 2;
        else
                return 0
        fi
}

GDAL_VERSION=($(gdal-config --version | tr '.' ' '))
echo "GDAL Version: `echo ${GDAL_VERSION[*]} | tr ' ' '.'`"

H_BAND=1
ALGORITHMS=()
while [[ $# > 0 ]];
do
        case $1 in
		--input)
			INPUT=$2
			shift
			shift
		;;
		--height)
			H_IMAGE=$2
			shift
			shift
		;;
		--output)
			OUTPUT=$2
			shift
			shift
		;;
		--height-band)
			H_BAND=$2
			shift
			shift
		;;
		--algorithms)
			eval "ALGORITHMS=($2)"
			shift
			shift
		;;
		*)
        	shift
		;;
        esac
done

if [ -z $INPUT ];
then
        echo [`date`] No input image specified >&2
        exit 2
fi

if [ -z $H_IMAGE ];
then
        echo [`date`] No height image specified >&2
        exit 2

fi

if [ -z $OUTPUT ];
then
        echo [`date`] No output file name specified >&2
        exit 2
fi


DEBUG_PREFIX="errcheck"
GEOS_HEIGHT_CORRECTION="${DEBUG_PREFIX} geosheightcorrection"
GDAL_2_XYZ="${DEBUG_PREFIX} gdal2xyz.py"
GDAL_TRANSLATE="${DEBUG_PREFIX} gdal_translate"
GDAL_EDIT="${DEBUG_PREFIX} gdal_edit.py"
GDAL_WARP="${DEBUG_PREFIX} gdalwarp"
GDAL_INFO="${DEBUG_PREFIX} gdalinfo"
GDAL_MERGE="${DEBUG_PREFIX} gdal_merge.py"
GDAL_GRID="${DEBUG_PREFIX} gdal_grid"
GDAL_CALC="${DEBUG_PREFIX} gdal_calc.py"
OGR_2_OGR="${DEBUG_PREFIX} ogr2ogr"
COPY="${DEBUG_PREFIX} gdalmanage copy"
MKDIR="${DEBUG_PREFIX} mkdir"
CD="${DEBUG_PREFIX} cd"
RM="${DEBUG_PREFIX} rm"
TEMP_DIR=/tmp/paralax_$$

$MKDIR $TEMP_DIR

BLANK_IMAGE=${TEMP_DIR}/blank.tif
REPR_H_IMAGE=${TEMP_DIR}/h_image.tif
$GDAL_CALC -A $INPUT --outfile=$BLANK_IMAGE --calc="0"
$GDAL_MERGE -separate -o $REPR_H_IMAGE $BLANK_IMAGE $BLANK_IMAGE
$GDAL_WARP -order 1 $H_IMAGE $REPR_H_IMAGE

COORDS_IMAGE=${TEMP_DIR}/coords_image.tif
$GEOS_HEIGHT_CORRECTION --input $REPR_H_IMAGE --height-band $H_BAND --output $COORDS_IMAGE

IMG_STAGE_1=${TEMP_DIR}/stage1.tif
${GDAL_MERGE} -separate -o ${IMG_STAGE_1} ${COORDS_IMAGE} $INPUT

BANDS_COUNT=`getBandsCount $INPUT`
IMG_STAGE_2=${TEMP_DIR}/stage2.csv
BAND_OPTIONS=`getXYZAllBandsOption $[$BANDS_COUNT + 2]`
${GDAL_2_XYZ} ${BAND_OPTIONS} ${IMG_STAGE_1} ${IMG_STAGE_2}

IMG_STAGE_2_1=${TEMP_DIR}/stage2_1.csv
BANDS_HEADER=`seq $[$BANDS_COUNT + 2] | sed 's#^#b#;s#$# #' | tr -d '\n' | sed 's# $##'` 
echo "x y ${BANDS_HEADER}" >${IMG_STAGE_2_1}
cut -f 3- -d' ' <${IMG_STAGE_2} >>${IMG_STAGE_2_1}

SRS_FILE=${TEMP_DIR}/srs.txt
getProjSRS ${IMG_STAGE_1} >${SRS_FILE}
IMG_STAGE_2_2=${TEMP_DIR}/stage2_2.shp
${OGR_2_OGR} -a_srs ${SRS_FILE}  -oo X_POSSIBLE_NAMES=x -oo Y_POSSIBLE_NAMES=y -oo AUTODETECT_TYPE=yes ${IMG_STAGE_2_2}  CSV:${IMG_STAGE_2_1} stage2_1


IMG_SIZE=`getSize $INPUT`
IMG_PIX_SIZE=`getPixelSize $INPUT`
IMG_PIX_SIZE_1=`echo $IMG_PIX_SIZE | cut -d " " -f 1`
IMG_PIX_SIZE_2=`echo $IMG_PIX_SIZE | cut -d " " -f 2`
IMG_CORNERS=`getULLR $INPUT`

IMG_STAGE_3=${TEMP_DIR}/stage3
for((B=1; B <= $BANDS_COUNT; B++))
do
	NO_DATA_V=`getNoData $INPUT $B`
	if [ "$NO_DATA_V" != "" ];
	then
		ALG_NODATA=":nodata=${NO_DATA_V}"
		NODATA_TEST="-where `getNoDataTest b${B} ${NO_DATA_V}`"	
	else
		ALG_NODATA=":nodata=nan"
	fi

	ALG=${ALGORITHMS[B - 1]}
	if [ -z $ALG ];
	then
		ALG=average
	fi

	${GDAL_GRID} -l stage2_2 -zfield b${B} $NODATA_TEST -outsize $IMG_SIZE -spat $IMG_CORNERS -a ${ALG}:radius1=${IMG_PIX_SIZE_1}:radius2=${IMG_PIX_SIZE_2}${ALG_NODATA} $IMG_STAGE_2_2 ${IMG_STAGE_3}_${B}.tif
done

IMG_STAGE_4=${TEMP_DIR}/stage4.tif
${GDAL_MERGE} -a_nodata nan -separate -o ${IMG_STAGE_4} ${IMG_STAGE_3}*.tif

$COPY $IMG_STAGE_4 $OUTPUT

# Alternative version of copying which preservs exact raster shape
#$COPY $INPUT $OUTPUT
#$GDAL_WARP -r near $IMG_STAGE_4 $OUTPUT

${RM} -r $TEMP_DIR
