"""
Methylation-based deconvolution of tissues-of-origin. Takes input files with tab-separated columns of methylation
levels (float between 0.0 and 1.0). Each line begins with chr - startpos - endpos - ID. Values start in column five.

One reference file and a variable number of sample files can be input, but coordinates MUST be identical for all
lines between all files
Copyright (C) 2020 by Florian Erger
Contact: florian.erger@uk-koeln.de
License: GNU GPLv3
"""

import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import minimize


def check_integrity(ref, samples):
    """
    Checks if reference and sample file contain the identical coordinates for deconvolution.
    """
    ref_hash = hash("".join(["".join(x[:3]) for x in ref]))
    for sample in samples:
        sample_hash = hash("".join(["".join(x[:3]) for x in sample]))
        if ref_hash != sample_hash:
            return False
    return True


def loss(x):
    """
    Loss function for least-squares minimization.
    Takes vector x as input.
    Returns squared difference between (dot product of x and the reference matrix) and the sample matrix
    """
    diff = []
    for bin in range(len(reference_matrix)):
        diff.append(np.dot(x, reference_matrix[bin]))
    diff = np.array(diff, dtype=float)
    return np.sum(np.square(diff - sample_matrix))


def deconvolute(ref_in, sample_in, sample_count, nans_in):
    if nans_in:
        if args.verbose: print("\nDetected {0} NaN fields in sample #{1}. Cleaning up input data... ".format(len(nans_in), sample_count+1), end="")

        for idx in nans_in:
            ref_in = np.delete(ref_in, idx, 0)
            sample_in = np.delete(sample_in, idx)

        if args.verbose: print("done!")

    global reference_matrix; reference_matrix = ref_in
    global sample_matrix; sample_matrix = sample_in
    res = minimize(loss, x0, method='SLSQP', constraints=cons,
                   bounds=[(0, 1) for j in range(ref_num)], options={'disp': args.verbose})
    deconv_results[0].append("'" + sample_names[sample_count] + "'")
    for idx, value in enumerate(res.x):
        deconv_results[idx+1].append(value)


parser = argparse.ArgumentParser(description="\nMethylation-based deconvolution of tissues-of-origin. Takes input files"
                                             " with tab-separated columns of methylation levels (float between 0.0 and "
                                             "1.0). Each line begins with chr - startpos - endpos - ID. Values start in"
                                             " column five. \n"
                                             "One reference file and a variable number of sample files can be input, "
                                             "but coordinates MUST be identical for all lines between all files!\n"
                                             "Files should not contain NaN.")
parser.add_argument("reference_file", type=argparse.FileType("r"),
                    help="File with methylation values for reference tissues.")
parser.add_argument("sample_files", type=argparse.FileType("r"), nargs="+",
                    help="One or more space separated paths to sample methylation files.")
parser.add_argument("--ref_header", choices=[True, False], type=bool,
                    help="Set to true if the reference file has a header. If not "
                         "specified, program will attempt to detect this automatically.")
parser.add_argument("--sample_header", choices=[True, False], type=bool,
                    help="Set to true if the sample file(s) has/have a header. If not "
                         "specified, program will attempt to detect this automatically.")
parser.add_argument("--verbose", action="store_true", help="Give this argument to display status update messages")
parser.add_argument("--ineq", action="store_true", help="Give this argument to use inequality constraint of 1 (sum of "
                                                        "components must be LESS THAN OR EQUAL to 1). If not specified,"
                                                        " equality constraint of 1 is used (sum of components must be "
                                                        "EQUAL to 1).")
parser.add_argument("--conversion_correction", action="store_true", help="Give this argument to attempt "
                                                                         "automatic correction for experimental "
                                                                         "over- or underconversion during the "
                                                                         "deconvolution calculations. Most useful "
                                                                         "without the --ineq option.")

args = parser.parse_args()

# OPEN REFERENCE BEDFILE AND LOAD INTO MEMORY
ref = [row for row in csv.reader(args.reference_file, delimiter="\t")]
args.reference_file.close()

if ("chr" not in ref[0][0]) and args.ref_header is None:
    if args.verbose: print("Header detected in reference file, setting --ref_header to true...\n")
    args.ref_header = True
elif args.ref_header is None:
    if args.verbose: print("No header detected in reference file, setting --ref_header to false...\n")
    args.ref_header = False

if args.ref_header:
    ref_organs = [x for x in ref[0][4:]]
    ref = ref[1:]
else:
    ref_organs = ["Tissue_ref_{0}".format(x) for x in range(1, int((len(ref[0]) - 3)))]
ref_num = len(ref_organs)

# OPEN PREPARED SAMPLE BEDFILE(S) AND LOAD INTO MEMORY
sample_names = []
samples = [[] for i in range(len(args.sample_files))]
for cnt, samplefile in enumerate(args.sample_files):
    samples[cnt] = [row for row in csv.reader(samplefile, delimiter="\t")]
    samplefile.close()

    if ("chr" not in ref[0][0]) and args.sample_header is None:
        if args.verbose: print("Header detected in sample file, setting --sample_header to true...\n")
        args.sample_header = True
    elif args.sample_header is None:
        if args.verbose: print("No header detected in sample file, setting --sample_header to false...\n")
        args.sample_header = False

    sample_names.append(Path(samplefile.name).name)
    if args.sample_header:
        samples[cnt] = samples[cnt][1:]

sample_num = len(sample_names)
if args.verbose: print("Files successfully loaded.\n")

# CHECK IF SAMPLE FILES HAVE SAME LOCI TEMPLATE AS REFERENCE FILE
if args.verbose: print("Checking coordinate identity... ", end="")
if not check_integrity(ref, samples):
    if args.verbose: print("\n\nReference methylation file and sample methylation coordinates not identical.\n"
                           "Make sure that all supplied files are divided into the identical .bed regions!\n\n"
                           "Aborting analysis...")
    exit(1)
if args.verbose: print("done!\n")

# ADD ARTIFICIAL REFERENCES FOR OVER- UNDERCONVERSION CORRECTION
if args.conversion_correction:
    ref_organs.append("Underconversion_correction")
    ref_organs.append("Overconversion_correction")
    ref_num += 2

    ref = [x + ["1.0"] + ["0.0"] for x in ref]

# DETECT REGIONS WITH UNKNOWN METHYLATION IN REFERENCE AND REMOVE THESE REGIONS ALSO FROM SAMPLE DATA
nans = [idx for idx, row in enumerate(ref) if "nan" in "_".join(row).lower()]
nans = sorted(nans, reverse=True)

if nans:
    if args.verbose: print("Detected {0} NaN fields in reference. Cleaning up... ".format(len(nans)), end="")

    for idx in nans:
        del ref[idx]
        for i in range(sample_num):
            del samples[i][idx]

    if args.verbose: print("done!\n")

# SET UP MINIMIZE
deconv_results = [["Sample"]] + [[x] for x in ref_organs]
cons = ({'type': 'ineq' if args.ineq else 'eq',
         'fun' : lambda x: 1.0 - np.sum(x)})
x0 = np.zeros(ref_num)
n_loci = []

for counter, sample in enumerate(samples):
    # DETECT REGIONS WITH UNKNOWN METHYLATION IN SAMPLE DATA AND REMOVE THESE REGIONS LATER
    nans = [idx for idx, row in enumerate(sample) if "nan" in "_".join(row).lower()]
    nans = sorted(nans, reverse=True)
    n_loci.append(len(ref) - len(nans))

    # CALCULATE TISSUES OF ORIGIN BASED ON SUPPLIED REFERENCE METHYLATIONS
    deconvolute(np.array([x[4:] for x in ref], dtype=float), np.array([x[4:] for x in sample], dtype=float).transpose(),
                counter, nans)

# OUTPUT RESULTS
print("\n\n")
for result in deconv_results:
    print("\t".join([str(x) for x in result]))

if sample_num > 1 and min(n_loci) != max(n_loci):
    loci_string = "between {0} and {1}".format(min(n_loci), max(n_loci))
else:
    loci_string = str(n_loci[0])

if args.verbose: print("\nTissue composition for all samples has been calculated!\n"
                       "A total of {0} loci were considered in the composition calculation.\n\n"
                       "Successfully finished... exiting".format(loci_string))
