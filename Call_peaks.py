"""
Call WPS peaks in a file containing normalized WPS values
Copyright (C) 2020 by Florian Erger
Contact: florian.erger@uk-koeln.de
License: GNU GPLv3
"""

import argparse
from pathlib import Path
import csv
import statistics as stats


parser = argparse.ArgumentParser(description="Automatically call WPS peaks in a normalized WPS file.")
parser.add_argument("input_file", help="File with normalized WPS values for peak calling.")

args = parser.parse_args()

WPS = csv.reader(open(Path(args.input_file), "r"), delimiter=" ")

VALCACHE = []
POSCACHE = []
PEAKPOS = []
CONTNEG = 0
for row in WPS:
    score = float(row[2])
    pos = ":".join(row[0:2])
    VALCACHE.append(score)
    POSCACHE.append(pos)

    if score < 0:
        CONTNEG += 1
    else:
        CONTNEG = 0

    if CONTNEG > 5:
        # get peak if 50<region<450, otherwise do nothing
        if 50 <= len(VALCACHE) <= 450:
            score_sum = 0
            sum_to_beat = 0
            pos_to_beat = []

            # calculate median of region and check which continuous window of values above the median has the
            # largest sum, take center of that window as peak

            median = stats.median(sorted(VALCACHE))
            for idx, val in enumerate(VALCACHE):
                if val >= median:
                    score_sum += val
                    PEAKPOS.append(POSCACHE[idx])
                else:
                    if score_sum > sum_to_beat:
                        sum_to_beat = score_sum
                        pos_to_beat = [x for x in PEAKPOS]
                        score_sum = 0
                        PEAKPOS = []

            print(pos_to_beat[len(pos_to_beat) // 2])

        # reset caches
        VALCACHE = []
        POSCACHE = []
        CONTNEG = 0

