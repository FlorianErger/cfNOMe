"""
Find closest distances between WPS peaks in two different WPS files
Copyright (C) 2020 by Florian Erger
Contact: florian.erger@uk-koeln.de
License: GNU GPLv3
"""

from pathlib import Path
from bisect import bisect_left
import argparse


def closest(comp_list, number):
    idx = bisect_left(comp_list, number)

    if idx == 0:
        return comp_list[0] - number
    if idx == len(comp_list):
        return comp_list[-1] - number

    return min(comp_list[idx] - number, comp_list[idx - 1] - number, key=abs)


parser = argparse.ArgumentParser(description="Find the distances (-500bp to 500bp) between the closest WPS peaks in two"
                                             " different WPS files.")
parser.add_argument("peak1", help="File with baseline WPS peak coordinates")
parser.add_argument("peak2", help="File with comparison WPS peak coordinates")

args = parser.parse_args()

PEAKS_1 = {}
PEAKS_2 = {}
CHROM = sorted(set([x.split(":")[0] for x in open(Path(args.peak2), "r").readlines()]))

occs = {}
for i in range(-500, 501):
    occs[i] = 0

for c in CHROM:
    PEAKS_1[c] = sorted([int(x.split(":")[1]) for x in open(Path(args.peak1), "r").readlines() if x.split(":")[0] == c])
    PEAKS_2[c] = sorted([int(x.split(":")[1]) for x in open(Path(args.peak2), "r").readlines() if x.split(":")[0] == c])

    for peak in PEAKS_1[c]:
        diff = closest(PEAKS_2[c], peak)
        if -500 <= diff <= 500:
            occs[diff] += 1

for i in range(-500, 501):
    print("{0}\t{1}".format(i, occs[i]))

