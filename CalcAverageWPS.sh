#!/bin/bash

# Average raw WPS values from file containing multiple blocks of values of identical length.
# Copyright (C) 2020 by Florian Erger
# Contact: florian.erger@uk-koeln.de
# License: GNU GPLv3

usage()
{
    echo -e "\nAverage raw WPS values from file containing multiple blocks of values of identical length.\n"
    echo -e "Usage: $0 [options] -i input_WPS_file --num number_of_blocks\n"
	echo -e "required arguments:\n\t-i FILEPATH\t\tWPS file containing blocks of scores of equal length in direct succession\n"
	echo -e "\t--num\t\t\tnumber of regions in the input file that need to be averaged"
    echo -e "optional arguments:\n\t--keep\t\t\tSet to 'true' to keep input file of raw WPS (default: false)\n"
	echo -e "\t-h, --help\t\tshow this help message and exit\n"
}

KEEP="false"

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
		--num)
			NUM_BLOCKS=${VALUE}
			;;
		--keep)
			KEEP=${VALUE}
			;;
        *)
            echo -e "\nERROR: unknown parameter \"${PARAM}\"\n"
            usage
            exit 1
            ;;
    esac
    shift 2
done

if [[ "${RAWFILE}" = "" || ${NUM_BLOCKS} = "" ]]
then
	echo -e "\nERROR: missing parameter(s)"
	usage
	exit 1
fi

scriptdir=$(dirname "$0")
BLOCK_LEN=$(( $(wc -l "${RAWFILE}" | cut -d' ' -f1) / ${NUM_BLOCKS} ))

echo -e "\nNumber of blocks specified: ${NUM_BLOCKS}, calculated block length: ${BLOCK_LEN}\n" > /dev/tty

echo "Creating normalized file..." > /dev/tty

# normalize by number of endpoints and ROIs, otherwise highly covered samples generate higher scores.
num_EP=$(awk '{sum_EP += $4}
	END {print sum_EP}' "${RAWFILE}")

awk -v EP=${num_EP} '
        {WPS[(NR-1)%'${BLOCK_LEN}'] += $3 / EP}
    END {for (i = 0; i < '${BLOCK_LEN}'; i++) {
            print i - int('${BLOCK_LEN}'/2), WPS[i]}
    }' "${RAWFILE}" > "${RAWFILE}"_avg_norm.txt

# additional file with running average for window of +-pmavg basepairs subtracted from consensus value
PMAVG=500
awk -v pmavg=${PMAVG} '
		{val[NR] = $2; sum += $2 - val[NR - 2*pmavg - 1]} NR > 2*pmavg {print $1 - pmavg, val[NR - pmavg] - (sum / (2*pmavg + 1))
	}' "${RAWFILE}"_avg_norm.txt > "${RAWFILE}"_avg_norm+-${PMAVG}_running_mean_zero.txt

if [[ ${KEEP} = "false" ]]
then
	rm "${RAWFILE}"
fi

