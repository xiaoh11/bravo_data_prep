#!/usr/bin/env python

import sys
import json
from statistics import mean, median

"""
Calculate statistics of aggregated pileups
Reads from stdin. Writes to stdout.
Expects tab delimited input in format chrom\tpos\tdepth;depth;depth...;depth
Final field, depth, is list of depths separated by ';'
E.g.
chr11   61400   11
chr11   61401   1;11
chr11   61402   3;2;1;5
chr11   61403   3;2;1;5
chr11   61404   3;1;2;1;5
chr11   61405   3;1;2;1;5
"""


def main(n_indiv=1000):
    depth_thresholds = [1, 5, 10, 15, 20, 25, 30, 50, 100]

    for line in sys.stdin:
        chrom, pos, depth_str = line.rstrip().split()
        depths = [int(i) for i in depth_str.split(';')]

        # Tabulate counts for each depth threshold to answer:
        #   How many depths are greater than or equal to each threshold?
        counts = [0] * len(depth_thresholds)
        for depth in depths:
            for idx, break_val in enumerate(depth_thresholds):
                if depth >= break_val:
                    counts[idx] += 1

        chromosome_id = chrom.replace('chr', '', 1)

        summary = {}
        summary["chrom"] = chromosome_id
        summary["start"] = pos
        summary["end"] = pos
        summary["mean"] = round(mean(depths), 4)
        summary["median"] = round(median(depths), 4)

        # Calculate proportion of individuals counted in each bin for current position
        for br, depth_count in zip(depth_thresholds, depths):
            proportion = depth_count / n_indiv
            summary[br] = proportion

        sys.stdout.write(f"{chromosome_id}\t{pos}\t{pos}\t{json.dumps(summary)}\n")


if __name__ == '__main__':
    main()
