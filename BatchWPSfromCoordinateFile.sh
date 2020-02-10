#!/bin/bash

# Calculate raw WPS values for genomic regions given in file and optionally average the values over all supplied regions.
# Copyright (C) 2020 by Florian Erger
# Contact: florian.erger@uk-koeln.de
# License: GNU GPLv3

usage()
{
    echo -e "\nCalculate WPS for multiple coordinate ranges given in a text file (one line per coordinate, e.g. 8:43546000-43548000).
Optionally, average the values over all supplied regions (for this, all regions must be the same length!)\n"
    echo -e "Usage: $0 [options] -i input_bam_file -c coordinate_file\n"
	echo -e "required arguments:\n\t-i FILEPATH\t\tsorted and indexed .bam file\n"
	echo -e "\t-c FILEPATH\t\ttext file with chromosomal coordinates to calculate WPS\n\n"
    echo -e "optional arguments:\n\t-w, --window\t\twindow length for WPS (default: 120)\n"
	echo -e "\t--min\t\t\tminimum length of fragments to be included (default: 1)\n"
	echo -e "\t--max\t\t\tmaximum length of fragments to be included (default: 500)\n"
	echo -e "\t--average\t\tset to 'false' if average over all regions should not be calculated and instead simple output of 
\t\t\t\traw WPS for all supplied regions should be generated (default: true)\n"
	echo -e "\t--keep\t\t\tset to 'true' if you want to keep the file with raw WPS values. Only relevant if --average
\t\t\t\tis active (default: false)\n"
	echo -e "\t-p\t\t\tNumber of parallel threads to be used. Job will be split by coordinate lines. You usually 
\t\t\t\twill want to do this! (default: 1)\n"
	echo -e "\t-h, --help\t\tshow this help message and exit\n"
}

WINDOW=120
MIN=1
MAX=500
AVERAGE="true"
KEEP="false"
CORES=1

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
            BAMFILE="${VALUE}"
            ;;
		-w|--window)
			WINDOW=${VALUE}
			;;
		-c)
			COORDFILE="${VALUE}"
			;;
		--min)
			MIN=${VALUE}
			;;
		--max)
			MAX=${VALUE}
			;;
		--average)
			AVERAGE=${VALUE}
			;;
		--keep)
			KEEP=${VALUE}
			;;
		-p)
			CORES=${VALUE}
			;;
        *)
            echo -e "\nERROR: unknown parameter \"${PARAM}\"\n"
            usage
            exit 1
            ;;
    esac
    shift 2
done

if [[ "${BAMFILE}" = "" || "${COORDFILE}" = "" ]]
then
	echo -e "\nERROR: missing parameter(s)"
	usage
	exit 1
fi

# if regions should be averaged, check if all regions are same size
if [[ ${AVERAGE} = "true" ]]
then
	echo -e "\nChecking coordinate file integrity...\n"
	while read line
	do
		len2=$(( $(echo ${line} | cut -d- -f2) - $(echo ${line} | cut -d- -f1 | cut -d: -f2) ))
		if [[ ${len1} != "" && ${len1} -ne ${len2} ]]
		then
			echo -e "\nERROR: --average is set to true but coordinate file contains regions of unequal length. Consider setting --average to false\n"
			exit 1
		fi
		len1=$(( $(echo ${line} | cut -d- -f2) - $(echo ${line} | cut -d- -f1 | cut -d: -f2) ))
	done < "${COORDFILE}"
fi

lines=$(wc -l "${COORDFILE}" | cut -d' ' -f1)

scriptdir=$(dirname "$0")
coordbase=$(basename "${COORDFILE}")

`parallel -k -j ${CORES} --eta "${scriptdir}"/CalcWPSforCoordinate.sh --window ${WINDOW} --min ${MIN} --max ${MAX} -i "${BAMFILE}" --region {} :::: "${COORDFILE}" >> "${BAMFILE}"_"${coordbase}"_"${WINDOW}"WPS_"${MIN}"-"${MAX}"bp`

if [[ ${AVERAGE} = "true" ]]
then
	`"${scriptdir}"/CalcAverageWPS.sh --keep ${KEEP} -i "${BAMFILE}"_"${coordbase}"_"${WINDOW}"WPS_"${MIN}"-"${MAX}"bp --num ${lines}`
fi

echo -e "\nCalculated average WPS for a window of ${WINDOW}bp for ${lines} coordinates."
echo -e "Only fragments between ${MIN}bp and ${MAX}bp were analysed."
echo -e "\nGenerated WPS output saved as ${BAMFILE}_${coordbase}_${WINDOW}WPS_${MIN}-${MAX}bp"

