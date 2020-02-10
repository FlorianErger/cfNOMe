The cfNOMe Toolkit

The following is an introduction into the bioinformatic tools used in the manuscript "cfNOMe – A single assay for comprehensive epigenetic analyses of cell-free DNA". These tools are meant to be used on paired-end
cytosine-converted (e.g. enzymatic cytosine conversion) cell-free DNA NGS data. In principle, these tools should be equally compatible with bisulfite-converted libraries.

1. CONTENTS

Included in this package are:

Software:

- CalcWPSforCoordinate.sh	bash script for calculating the WPS given an alignment file and a chromosomal location.
- BatchWPSfromCoordinateFile.sh	bash script for calculating the WPS in batch given an alignment file and a file containing any number of chromosomal locations (in shape 8:43540000-43550000).
- CalcAverageWPS.sh		bash script for averaging and normalizing raw WPS data, especially when trying to merge signals from multiple locations (e.g. transcription start sites of multiple genes).
- Call_peaks.sh			bash script for calling peaks - representing nucleosome or transcription factor binding positions - from raw WPS data.
- Call_peaks.py			Python3 program being automatically called by Call_peaks.sh for performing the peak calling logic.
- find_nearest_peak.py		Python3 program that takes two input files containing chromosomal peak locations and outputs the distance between peak calls in the two samples. Distances can be -500 to +500bp.
- methylation_deconvolution.py	Python3 program that takes tab-separated methylation data as input (reference methylation and sample methylation are both needed) and calculates the most likely sample composition.

2. DEPENDENCIES

GNU parallel v20161222
samtools v1.7
Python v3.6.8
 - numpy v1.16.1
 - scipy v1.2.0
bismark v0.20.0
bedmap v2.4.26
sort-bed v2.4.26
bowtie2 v2.3.4.1

3. INSTALLATION

No special installation should be necessary. Just clone the repo onto your harddrive:

git clone https://github.com/FlorianErger/cfNOMe 

and execute the scripts/programs from the command line. Usage instructions for each component of the toolkit are available using the --help command line parameter.

4. PREPARING INPUT DATA FOR YOUR OWN ANALYSES

To perform cfNOMe analyses on your own data, perform the following steps first:

a. Perform methylation sensitive alignment using bismark (e.g. bismark --multicore * --output_dir /path/to/output/directory --genome /path/to/bismark/prepared/reference -1 /path/to/read1.fq.gz -2 /path/to/read2.fq.gz)
b. Perform methylation summary analysis on output .bam file from step a (e.g. bismark_methylation_extractor --parallel * --comprehensive --output /path/to/output/directory --merge_non_CpG --bedGraph --gzip /path/to/bismark.bam)
c. sort and index the bismark output .bam file from step a (e.g. samtools sort -@ * /path/to/input.bam > input_sorted.bam && samtools index -@ * input_sorted.bam)
d. Generate methylation summaries for loci of interest from bismark methylation data from step b (e.g. gunzip -c /path/to/bismark.bedGraph.gz | awk '{if(NR>1){print "chr"$1"\t"$2"\t"$3"\tid-"NR-1"\t"$4/100}}' | sort-bed - | bedmap --echo --delim '\t' --mean /path/to/CpG_loci_to_be_studied.bed - > average_methylation_at_loci_of_interest.bed)
*: enter number of parallel threads to be used

For a cfDNA dataset sequenced at ~30X, these steps should take about 8h for a, 1h for b, 1h for c, 1h for d on a normal desktop PC, provided appropriate parallelization settings are used.

5. PERFORMING cfNOMe ANALYSES (DEMO DATA)

Included with this toolkit are demo datasets. For the Nucleosome footprinting tools, a sorted and indexed .bam file has been included, containing NGS reads on hg19 chr8:43540000-43550000. To perform analyses on this
demo data, use the following commands:

a. To calculate raw WPS in an area of interest: 											./CalcWPSforCoordinate.sh -i Example_data/Nucleosome_footprinting/NF_demo_data.bam -r 8:43546000-43548000 > Raw_WPS.txt
b. Optional: normalize the raw WPS (incurs information loss at both ends of the window):	./CalcAverageWPS.sh --keep true -i Raw_WPS.txt --num 1
c. Optional: call WPS peaks:																./Call_peaks.sh -i Raw_WPS.txt > Peak_calls.txt
d. Optional: compare peak calls (requires second sample):									python3 find_nearest_peak.py Peak_calls_from_sample_1.txt Peak_calls_from_sample_2.txt

The expected run time for each of these steps is <1 second on a normal desktop PC. For bigger analyses with larger alignment files and larger regions of interest (over several Mb), step a can be parallelized with the -p flag. Runtime for a 200Mb sequence chunk is estimated at ~1 hour for steps a-d.

For methylation deconvolution analysis, we include the reference methylation levels at 7890 loci-of-interest from the publication by Moss et al., Nat Commun., 2018 for 25 reference tissues and for one of our healthy individuals. To 
perform analyses on this demo data, use the following commands:

e. calculate tissue contribution based on methylation levels:								python3 methylation_deconvolution.py Example_data/Methylation/Reference_tissue_methylation.bed Example_data/Methylation/Sample_methylation_levels.bed --ineq > Tissue_fractions.tsv

The expected run time for a methylation deconvolution using these references and loci-of-interest is on the order of 10 seconds on a normal desktop PC.

Expected output files can be found in the respective subfolders of the Example_data directory.

6. PERFORMING BATCH WPS CALCULATIONS OVER MULTIPLE COORDINATES

For calculating the WPS at several specific coordinates (and optionally averaging them to a single value, e.g. for looking at the transcription start sites of several thousand genes at once), use:

a. Use coordinate file to calculate WPS over many different regions and output averaged result:	./BatchWPSfromCoordinateFile.sh --average true --keep true -i /path/to/your_sorted_and_indexed.bam -c Example_data/Nucleosome_footprinting/NF_TSS_coordinates_for_10k_genes+-5kb.txt
