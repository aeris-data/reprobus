#!/bin/bash

set -e

SCRIPT_NAME=$(basename "$0")

function help() {
    bold=$(tput bold)
    normal=$(tput sgr0)
    echo ""
    echo "###   SCRIPT FOR            _(  )_( )_    "
    echo "###    THE ECMWF           (_   _    _)   "
    echo "###     DATA EXTRACTION   / /(_) (__)     "
    echo "###                      / / / / / /      "
    echo "###                     / / / / / /       "
    echo "###"
    echo "### ${bold}${SCRIPT_NAME}${normal} script handles the ECMWF O3 data extraction needed for the"
    echo "### REPROBUS simulations. The user have to define the start and end date of the data"
    echo "### extraction and the paths to the working directory and data directory."
    echo "### The working directory will just contain the request file for the MARS API; the data"
    echo "### directory will contain the final extracted data; the working and data directories"
    echo "### can be the same. The data will be extracted on the model levels of the ECMWF."
    echo "### The extraction is performed on the MARS server, hence the script must be launched"
    echo "### on MARS server."
    echo "###"
    echo "### Usage: ${SCRIPT_NAME} [options] arguments"
    echo "### Options:"
    echo "###   ${bold}-h, --help${normal}              Show this help message and exit"
    echo "### Arguments:"
    echo "###   ${bold}--config conf_filepath${normal}  This argument must correspond to the configuration"
    echo "###                           file where the user defines input parameters needed"
    echo "###                           for the extraction."
    echo "###"
    echo "### +--------------------------------------------------------------------------------+"
    echo "### | Example filename : ${bold}extraction.conf${normal}                                             |"
    echo "### | Example content below :                                                        |"
    echo "### |                                                                                |"
    echo "### | START_DATE='20230101'                                                          |"
    echo "### | DATA_DIR='/my/dir/for/data'                                                    |"
    echo "### | WORKING_DIR='/working/dir/for/aux/files'                                       |"
    echo "### +--------------------------------------------------------------------------------+"
    echo ""
}

function info_msg(){
    txt=$1
    echo "$(date +'%d/%m/%Y %H:%M:%S')   [INFO]   ${txt}"
}

function err_msg(){
    txt=$1
    echo "$(date +'%d/%m/%Y %H:%M:%S')   [ERROR]   ${txt}"
}

function warn_msg(){
    txt=$1
    echo "$(date +'%d/%m/%Y %H:%M:%S')   [WARNING]   ${txt}"
}

function check_args(){
    if [[ -z ${START_DATE} ]]; then
        err_msg "No start date was defined. Exiting script."
        exit 1
    fi
    if [[ -z ${DATA_DIR} ]]; then
        err_msg "No data directory was defined. Exiting script."
        exit 1
    fi
    if [[ -z ${WORKING_DIR} ]]; then
        err_msg "No working directory was defined. Exiting script."
        exit 1
    fi
    if [[ ${START_DATE} =~ ^[0-9]{8}$ ]]; then
        exit_status=0
    else
        err_msg "Starting date 'START_DATE' is not in the correct format YYYYMMDD"
        exit_status=1
    fi
    return ${exit_status}
}

function write_fortran_script(){
    # nb_loop=138 in this Fortran script because we have :
    #  - 1 variable on 137 levels (ozone) + 1 variable on single level (lnsp)
    #  - each variable is extracted on times 0,3,6,9,12,15,18,21 hour of the day (8 different times)
    #  - (1x137 + 1)x8 = 5488 records for one day of extracted data
    # and each record is on lat-lon grid, so nb_data=(nx)x(ny) pixels
    nx=$(bc <<< "scale=0; 360/${SPATIAL_RESOLUTION}")
    ny=$(bc <<< "scale=0; 180/${SPATIAL_RESOLUTION} + 1")
    ndata=$((${nx}*${ny}))
    cat > ${WORKING_DIR}/grib_rep.f90 << EOF
PROGRAM grib_rep
    implicit none
    CHARACTER(LEN=*), PARAMETER        :: input_file='data'
    CHARACTER(LEN=*), PARAMETER        :: header_file='entete'
    CHARACTER(LEN=26)                  :: header
    INTEGER, PARAMETER                 :: nb_loop=138
    INTEGER, PARAMETER                 :: nb_data=${ndata}
    integer                            :: err1,err2
    integer                            :: year,month,day,hour,stp1,c
    integer                            :: vers,par,ilev,lev,ctre,minu
    integer                            :: gpi,gdef,utr,stp2,tri,flg,i,j
    real*8                             :: lat,lon,tab(nb_data)
    !**********************************************!
    !    Open header and data file for reading     !
    !**********************************************!
    flg=128
    c=0
    ilev=109
    ctre=98
    open(11,file=trim(header_file),iostat=err1,status='old')
    if (err1 /= 0) print*, "Attention le fichier header n existe pas."
    open(25,file=trim(input_file),iostat=err2,status='old')
    if (err2 /= 0) print*, "Attention le fichier d entree n existe pas."

    do i=1,nb_loop
        read(11,*,iostat=err1) gpi,gdef,par,lev,year,month,day,hour,minu,utr,stp1
        write(10) flg,ctre,gpi,gdef,flg,par,ilev,lev,c,year,month,day,hour,minu,utr,stp1,c,c
        year=year-2000
        read(25,*,iostat=err2) header
        do j=1,nb_data
            read(25,*,iostat=err2) lat,lon,tab(j)
        enddo
        write(10)(tab(j),j=1,nb_data)
    enddo
    close(11)
    close(25)
    close(10)
END PROGRAM grib_rep
EOF
    gfortran -march=native -O3 -Wno-tabs -fconvert=big-endian -o ${WORKING_DIR}/grib_rep ${WORKING_DIR}/grib_rep.f90
}

function main(){
    module load ecmwf-toolbox

    SPATIAL_RESOLUTION="2."

    write_fortran_script
    cp ${WORKING_DIR}/grib_rep ${DATA_DIR}/

    DATE=${START_DATE}

    cat > ${WORKING_DIR}/data.req <<EOF
retrieve,
    class    = od,
    stream   = oper,
    date     = ${DATE},
    type     = an,
    time     = 12,
    levtype  = ml,
    levels   = 1,
    area     = 90./0./-90./360.,
    resol    = 106,
    grid     = ${SPATIAL_RESOLUTION}/${SPATIAL_RESOLUTION},
    param    = lnsp,
    target   = "${DATA_DIR}/datafile"
retrieve,
    levels   = 1/to/137/by/1,
    param    = o3,
    target   = "${DATA_DIR}/datafile"
EOF

    mars ${WORKING_DIR}/data.req

    grib_get -p generatingProcessIdentifier,bitMapIndicator,paramIdECMF,level,year,month,day,hour,minute,indicatorOfUnitOfTimeRange,forecastTime ${DATA_DIR}/datafile > ${DATA_DIR}/entete
    grib_get_data -M -F "%.15e" ${DATA_DIR}/datafile > ${DATA_DIR}/data

    cd ${DATA_DIR}
    ${DATA_DIR}/grib_rep

    mv ${DATA_DIR}/fort.10 ${DATA_DIR}/ecmwf_o3_${DATE}

    rm ${DATA_DIR}/data ${DATA_DIR}/entete ${DATA_DIR}/datafile

}

# ----------------------------------------------------------------------------------------------
# BASH SCRIPT
# ----------------------------------------------------------------------------------------------

opts=$(getopt --longoptions "config:,help" --name "$(basename "$0")" --options "h" -- "$@")
eval set --$opts

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) shift; CONFIG_FILE=$1; shift;;
        -h|--help) help; exit 0; shift;;
        \?) shift; err_msg "Unrecognized options"; exit 1; shift;;
        --) break;;
    esac
done

if [[ -z ${CONFIG_FILE} ]]; then
    err_msg "No configuration file was passed. Exiting script."
    exit 1
fi

source ${CONFIG_FILE}
check_args
ARGS_STATUS=$?
if [[ ${ARGS_STATUS} == 1 ]]; then
    err_msg "Errors in input parameters detected, job was not launched"
    exit 1
else
    main
fi