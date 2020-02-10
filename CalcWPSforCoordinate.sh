#!/bin/bash

# Calculate raw WPS values in a given genomic range.
# Copyright (C) 2020 by Florian Erger
# Contact: florian.erger@uk-koeln.de
# License: GNU GPLv3

usage()
{
    echo -e "\nCalculate WPS for a given coordinate range\n"
    echo -e "Usage: $0 [options] -i input_bam_file -r region\n"
	echo -e "required arguments:\n\t-i FILEPATH\t\tsorted and indexed .bam file\n"
	echo -e "\t-r, --region\t\tchromosomal coordinates to calculate WPS in (e.g. 8:43546000-43548000)\n\n"
    echo -e "optional arguments:\n\t-w, --window\t\twindow length for WPS (default: 120)\n"
	echo -e "\t--min\t\t\tminimum length of fragments to be included (default: 1)\n"
	echo -e "\t--max\t\t\tmaximum length of fragments to be included (default: 500)\n"
	echo -e "\t-p\t\t\tNumber of parallel threads to be used. Job will be split into 1Mb chunks. You usually will want to 
\t\t\t\tdo this! (default: 1)\n"
	echo -e "\t-h, --help\t\tshow this help message and exit\n"
}

WINDOW=120
MIN=1
MAX=500
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
		-r|--region)
			CHR=$(echo ${VALUE} | cut -d: -f1)
			STARTPOS=$(echo ${VALUE} | cut -d: -f2 | cut -d- -f1)
			ENDPOS=$(echo ${VALUE} | cut -d: -f2 | cut -d- -f2 | cut -d$'\r' -f1)
			;;
		--min)
			MIN=${VALUE}
			;;
		--max)
			MAX=${VALUE}
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

if [[ "${BAMFILE}" = "" || ${CHR} = "" ]]
then
	echo -e "\nERROR: missing parameter(s)"
	usage
	exit 1
fi

scriptdir=$(dirname "$0")

if [[ $(( ENDPOS - STARTPOS )) -gt 1000000 ]]
then
	BINS=$(( (ENDPOS-STARTPOS) / 1000000 ))
	printf "\n\nSplitting job into "$(( ${BINS} + 1 ))"x1Mb bins for parallelisation...\n\n" > /dev/tty
	STARTARR=($(seq 0 ${BINS} | awk -v OFS=" " '
			{if ('${STARTPOS}' + (($1 + 1) * 1000000) < '${ENDPOS}') {
				print '${STARTPOS}' + ($1 * 1000000) "-" '${STARTPOS}' + (($1 + 1) * 1000000)}
			else {
				print '${STARTPOS}' + ($1 * 1000000) "-" '${ENDPOS}'}
			}'))
	`parallel -k -j ${CORES} --eta "${scriptdir}"/CalcWPSforCoordinate.sh --window ${WINDOW} --min ${MIN} --max ${MAX} -i "${BAMFILE}" --region ${CHR}:{} ::: ${STARTARR[@]} >> "${BAMFILE}"_TEMP`
	cat "${BAMFILE}"_TEMP; rm "${BAMFILE}"_TEMP
	exit
fi

#           -F 0x10 flag extracts only forward reads, because we don't care about the read ends in PE, but the
#           fragment ends! We then add the fragment length to get the fragment end point. We extract start-500
#			to make sure that we have no edge artefacts when calculating the WPS. Scores should be for window
#			centered on position in question, so adjustments with half window lengths have to be made.

saminput=$(samtools view -F 0x10 "${BAMFILE}" "${CHR}":"$((STARTPOS-500))"-"$((ENDPOS+200))" | \
		   awk '$9>='${MIN}' && $9<='${MAX}'')
TEMPORARY+=$(echo "${saminput}" | awk -v chr=${CHR} '
BEGIN   {OFS="\t"; WindowLength = '${WINDOW}'}

        {Endpoints[$4]++; Endpoints[$4 + $9 - 1]++}
        {if (NR == 1) {FirstPos = $4}}
        {for (WindowCenter=$4 + (WindowLength/2) + 1; WindowCenter<$4 + $9 - (WindowLength/2); WindowCenter++) {
            ReadsStretchingWindow[WindowCenter]++}
        }

END     {LastPos = $4; if ('${STARTPOS}' != 0) {FirstPos = '${STARTPOS}'; LastPos = '${ENDPOS}'} 
        for (CurrentPos=FirstPos; CurrentPos<LastPos; CurrentPos++) {
            for (j=-WindowLength/2; j<WindowLength/2; j++) {
                if (CurrentPos+j in Endpoints) {EndpointsInWindow[CurrentPos] = \
                                                EndpointsInWindow[CurrentPos] + Endpoints[CurrentPos + j]}
                }
            print "chr"chr, CurrentPos, ReadsStretchingWindow[CurrentPos] - EndpointsInWindow[CurrentPos], Endpoints[CurrentPos]!="" ? Endpoints[CurrentPos] : "0"}
        }
')

echo "${TEMPORARY}"
