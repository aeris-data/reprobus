#!/bin/bash

set -e
SCRIPT_PATH=$0
SCRIPT_NAME=$(basename ${SCRIPT_PATH})

function help() {
    bold=$(tput bold)
    normal=$(tput sgr0)
    echo ""
    echo "###       .-'';'-."
    echo "###     ,'   <_,-.\`.      REPROBUS tool for"
    echo "###    /)   ,--,_>\_\       reactive processes"
    echo "###   |'   (      \_ |       ruling the ozone budget"
    echo "###   |_    \`-.    / |         in the stratosphere"
    echo "###    \\\`-.   ;  _(\`/"
    echo "###     \`.(    \/ ,' "
    echo "###       \`-....-'"
    echo "###"
    echo "### This script handles the REPROBUS simulation input"
    echo "### parameters and launches the simulation, as well as"
    echo "### the post-processing for the output additional"
    echo "### formatting and results visualization"
    echo "###"
    echo "### Usage: ${SCRIPT_NAME} [options] arguments"
    echo "### Options:"
    echo "###   ${bold}-h, --help${normal}                 Show this help message and exit"
    echo "### Arguments:"
    echo "###   ${bold}--config conf_filepath${normal}     This argument must correspond to the configuration"
    echo "###                              file where the user defines input parameters needed"
    echo "###                              for the simulation"
    echo "###"
    echo "### +--------------------------------------------------------------------------------+"
    echo "### | Example filename : ${bold}simulation.conf${normal}                                             |"
    echo "### | Example content below :                                                        |"
    echo "### |                                                                                |"
    echo "### | EXP=\"001442\"                                                                   |"
    echo "### | COMPILER=\"ifort\"                                                               |"
    echo "### | WDIR=\"/sedoo/resos/reprobus/WORK\"                                              |"
    echo "### | DATA_DIR=\"/sedoo/resos/reprobus/ECMWF\"                                         |"
    echo "### | RESTART_DIR=\"/sedoo/resos/reprobus/RESTART\"                                    |"
    echo "### | RES_DIR=\"/sedoo/resos/reprobus/RES\"                                            |"
    echo "### | IMAGES_DIR=\"/sedoo/resos/reprobus/IMAGES\"                                      |"
    echo "### | START_DATE=\"20230922\"                                                          |"
    echo "### | END_DATE=\"20230923\"                                                            |"
    echo "### | NSTART=1                                                                       |"
    echo "### +--------------------------------------------------------------------------------+"
    echo ""
}

function info(){
	txt=$1
	echo "$(date +'%d/%m/%Y %H:%M:%S')   [INFO]   ${txt}"
}
function error(){
	txt=$1
	echo "$(date +'%d/%m/%Y %H:%M:%S')   [ERROR]   ${txt}"
}
function warning(){
	txt=$1
	echo "$(date +'%d/%m/%Y %H:%M:%S')   [WARNING]   ${txt}"
}

function check_args(){
    if [[ "${COMPILER}" != "ifort" && "${COMPILER}" != "nvidia" ]]; then
        error "Unrecognized compiler options. Options are 'nvidia' or 'ifort', check your configuration file"
        exit_status=1
    fi
    if [[ ${START_DATE} =~ ^[0-9]{8}$ ]]; then
        exit_status=0
    else
        error "Starting date is not in the correct format YYYYMMDD"
        exit_status=1
    fi
    if [[ ${END_DATE} =~ ^[0-9]{8}$ ]]; then
        exit_status=0
    else
        error "Ending date is not in the correct format YYYYMMDD"
        exit_status=1
    fi
    return ${exit_status}
}

function count_days(){
    start_timestamp=$(date -d "${START_DATE}" +%s)
    end_timestamp=$(date -d "${END_DATE}" +%s)
    time_diff=$((end_timestamp - start_timestamp))
    days_diff=$((time_diff / 86400))
    echo ${days_diff}
}

function cleanup(){
    if [ ! -z ${WDIR} ] && [ -d ${WDIR} ]; then
        rm -f ${WDIR}/*
    fi
}

function prepare_directories(){
    if [ ! -d ${WDIR} ]; then mkdir -p ${WDIR}; fi
    if [ ! -d ${RESTART_DIR} ]; then mkdir -p ${RESTART_DIR}; fi
    if [ ! -d ${RES_DIR} ]; then mkdir -p ${RES_DIR}; fi
    if [ ! -d ${IMAGES_DIR} ]; then mkdir -p ${IMAGES_DIR}; fi
}

function main(){

    info "!======================================================================!"
    info "!        _____  ______ _____  _____   ____  ____  _    _  _____        !"
    info "!       |  __ \|  ____|  __ \|  __ \ / __ \|  _ \| |  | |/ ____|       !"
    info "!       | |__) | |__  | |__) | |__) | |  | | |_) | |  | | (___         !"
    info "!       |  _  /|  __| |  ___/|  _  /| |  | |  _ <| |  | |\___ \        !"
    info "!       | | \ \| |____| |    | | \ \| |__| | |_) | |__| |____) |       !"
    info "!       |_|  \_\______|_|    |_|  \_\\\____/|____/ \____/|_____/        !"
    info "!======================================================================!"

    prepare_directories
    cleanup
    cd ${WDIR}

    info "Compiling REPROBUS"

    ln -s /usr/local/REPROBUS/src/reprobus_l137_20230723.f ${WDIR}/reprobus_${EXP:2:4}.f
    ln -s /usr/local/REPROBUS/src/jno.f90 ${WDIR}/jno.f90
    ln -s /usr/local/REPROBUS/src/altitude.f90 ${WDIR}/altitude.f90
    ln -s /usr/local/REPROBUS/src/h2so4_l137_195101-200512.20140212 ${WDIR}/h2so4.txt
    ln -s /usr/local/REPROBUS/src/mopitt_corrected.txt ${WDIR}/mopitt_corrected.txt
    ln -s /usr/local/REPROBUS/src/jstrato.20111013 ${WDIR}/jstrato.txt
    ln -s /usr/local/REPROBUS/src/relief.txt ${WDIR}/relief.txt
    ln -s /usr/local/REPROBUS/src/ecmwf_137_levels.txt ${WDIR}/ecmwf_137_levels.txt

    if [[ "${COMPILER}" == "ifort" ]]; then
        source /opt/intel/oneapi/setvars.sh
        ifx -r8 -mcmodel=medium -convert big_endian -o reprobus_1442 reprobus_1442.f jno.f90 altitude.f90
    fi
    if [[ "${COMPILER}" == "nvidia" ]]; then
        /opt/nvidia/hpc_sdk/Linux_x86_64/23.9/compilers/bin/pgf90 -C -r8 -mcmodel=medium -fast -Kieee -byteswapio -o reprobus_${EXP:2:4} reprobus_${EXP:2:4}.f jno.f90 altitude.f90
    fi

    N_DAYS=$(count_days)
    DATE_IN=${START_DATE}
    DATE_OUT=${END_DATE}
    info "Simulation will be performed through ${NDAYS} days, from ${START_DATE} to ${END_DATE}"

    while [ ${DATE_IN} -lt ${DATE_OUT} ]; do

        date_start=${DATE_IN}
        date_end=$(date -d "${DATE_IN} + 1 day" "+%Y%m%d")

        info "Preparing input data for the simulation"

        ln -s ${DATA_DIR}/ecmwf_${date_start} ${WDIR}/ecmwf_${date_start}
        ln -s ${DATA_DIR}/ecmwf_${date_end} ${WDIR}/ecmwf_${date_end}

        if [[ ${NSTART} -eq 1 ]]; then
            ln -s ${RESTART_DIR}/MODEL_history_${date_start}12_${EXP} ${WDIR}/fort.90
        else
            # qinit_file=$(find /usr/local/REPROBUS/src -iname "qinit2d*")
            # ln -s ${qinit_file} ${WDIR}/qinit2d.txt
            # ln -s ${DATA_DIR}/ecmwf_o3_20240501 ${WDIR}/ecmwf_o3_20240501
            warning "This case is not yet implemented, you might want to contact the model providers for more information"
            warning "Exiting the script"
            exit 1
        fi

        info "Working directory for the current date simulation"
        ls -l ${WDIR}

        info "Launching REPROBUS"

        ./reprobus_${EXP:2:4}
        status=$?

        if [[ ${status} == 0 ]]; then
            info "Simulation successful, continuing with the script"
        else
            error "Something went wrong with the simulation, check the log messages above for more information"
            exit 1
        fi

        info "Organizing output files"

        mv history ${RESTART_DIR}/MODEL_history_${date_end}12_${EXP}
        mv stations.txt ${RES_DIR}/stations_${date_end}12_${EXP}
        mv fort.60 ${RES_DIR}/reprobus_kiruna_${date_end}12_${EXP}
        mv fort.61 ${RES_DIR}/reprobus_ohp_${date_end}12_${EXP}
        mv fort.62 ${RES_DIR}/reprobus_nyaalesund_${date_end}12_${EXP}
        mv fort.63 ${RES_DIR}/reprobus_sodankyla_${date_end}12_${EXP}
        mv fort.64 ${RES_DIR}/reprobus_yakutsk_${date_end}12_${EXP}
        mv fort.65 ${RES_DIR}/reprobus_ddu_${date_end}12_${EXP}
        mv fort.66 ${RES_DIR}/reprobus_harestua_${date_end}12_${EXP}
        mv fort.67 ${RES_DIR}/reprobus_marambio_${date_end}12_${EXP}
        mv fort.68 ${RES_DIR}/reprobus_southpole_${date_end}12_${EXP}
        mv fort.69 ${RES_DIR}/reprobus_airesadour_${date_end}12_${EXP}
        mv fort.70 ${RES_DIR}/reprobus_eureka_${date_end}12_${EXP}
        mv fort.71 ${RES_DIR}/reprobus_niamey_${date_end}12_${EXP}
        mv fort.72 ${RES_DIR}/reprobus_teresina_${date_end}12_${EXP}

        python3 /home/resos/git/reprobus/post-process-reprobus.py --date ${date_end} --restart-dir ${RESTART_DIR} --res-dir ${RES_DIR} --image-dir ${IMAGES_DIR}

#        rm ${WDIR}/ecmwf_${date_start}
#        rm ${WDIR}/ecmwf_${date_end}
#        if [[ ${NSTART} == 1 ]]; then
#            rm ${WDIR}/fort.90
#        else
#            rm ${WDIR}/qinit2d.txt
#            rm ${WDIR}/ecmwf_o3_20221201
#        fi

        DATE_IN=$(date -d "${DATE_IN} + 1 day" "+%Y%m%d")

    done
}

# +----------------------------------+
# | BASH SCRIPT                      |
# +----------------------------------+

opts=$(getopt --longoptions "config:,help" --name "$(basename "$0")" --options "h" -- "$@")
eval set --$opts

while [[ $# -gt 0 ]]; do
	case "$1" in
		--config) shift; CONFIG_FILE=$1; shift;;
        -h|--help) help; exit 0; shift;;
		\?) shift; error "Unrecognized options"; exit 1; shift;;
		--) break;;
	esac
done

if [[ -z ${CONFIG_FILE} ]]; then
    error "No configuration file was passed. Exiting script..."
    exit 1
fi

source ${CONFIG_FILE}
check_args
ARGS_STATUS=$?
if [[ ${ARGS_STATUS} == 1 ]]; then
    error "Errors in input parameters detected, job was not launched"
    exit 1
else
    main
fi
