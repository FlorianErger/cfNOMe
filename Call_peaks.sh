#!/bin/bash

# Calculate normalized WPS values and call Call_peaks.py to make WPS peak calls.
# Copyright (C) 2020 by Florian Erger
# Contact: florian.erger@uk-koeln.de
# License: GNU GPLv3

usage()
{
    echo -e "\nMake WPS peak calls from raw WPS files. Requires Call_peaks.py\n"
    echo -e "Usage: $0 [options] -i input_WPS_file\n"
	echo -e "required arguments:\n\t-i FILEPATH\t\tRaw WPS file for peak calling\n\n"
    echo -e "optional arguments:\n\t-h, --help\t\tshow this help message and exit\n"
}

while [[ "$1" != "" ]]
do
    PARAM=$(echo $1)
    VALUE=$(echo $2)
    case ${PARAM} in
        -h|--help)
            usage
            exit
            ;;
        -i)
            RAWFILE="${VALUE}"
            ;;
        *)
            echo -e "\nERROR: unknown parameter \"$PARAM\"\n"
            usage
            exit 1
            ;;
    esac
    shift 2
done

if [[ "${RAWFILE}" = "" ]]
then
	echo -e "\nERROR: missing parameter(s)"
	usage
	exit 1
fi

scriptdir=$(dirname "$0")

awk -v pmavg=500 '
	{val[NR] = $3; sum += $3 - val[NR - 2*pmavg - 1]} NR > 2*pmavg {print $1, $2 - pmavg, val[NR - pmavg] - (sum / (2*pmavg + 1))
}' "${RAWFILE}" > "${RAWFILE}"_rngavg_1000bp.txt

echo -e "\nTemporary normalized WPS file created, outputting peak calls...\n" > /dev/tty

echo $(python3 "${scriptdir}"/Call_peaks.py "${RAWFILE}"_rngavg_1000bp.txt) | tr ' ' '\n'

echo -e "\nPeak calling complete. Deleting temporary files...\n" > /dev/tty

rm "${RAWFILE}"_rngavg_1000bp.txt
