#!/usr/bin/env python

# This version of the depth statistic binning script handles an upstream issue where the positions
#  where not serialized as integers in the json data.  The upstream correction has been made.
#  This script exists to handle the interrim creation of testing data set while prod data runs.

import argparse
import gzip
import pysam
import rapidjson
import sys

argparser = argparse.ArgumentParser(description = 'Prunes base coverage by grouping bases into bins of similar median and mean depths. Starting with base X that has mean and median depths Z(X) and Y(X), every next base X+1 is added to the same bin if |Z(X+1) - Z(X)| <= LIMIT and |Y(X+1) - Y(X)| <= LIMIT')
argparser.add_argument('-i', '--in', metavar = 'file', dest = 'in_coverage_file', required = True, help = 'Input JSON coverage (compressed with gzip/bgzip) file')
argparser.add_argument('-l', '--limit', metavar = 'float', dest = 'fluctuation_limit', required = True, type = float, help = 'Threshold for maximal fluctuation of median and mean depths within a bin.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_coverage_file', required = True, help = 'Output JSON coverage (compressed with bgzip) file')


def write(output_file, data):
    # Cast to integer to handle upstream error serializing positions as strings.
    output_file.write('{}\t{:d}\t{:d}\t'.format(data['chrom'], data['start'], data['end']).encode())
    rapidjson.dump(data, output_file)
    output_file.write('\n'.encode())


if __name__ == '__main__':
    args = argparser.parse_args()
    with gzip.open(args.in_coverage_file, 'rt') as ifile, pysam.BGZFile(args.out_coverage_file, 'w') as ofile:
        line = ifile.readline()
        if not line:
            sys.exit(0)
        chrom, start, stop, data = line.rstrip().split('\t')

        # Use current line's data as bin initially and after each bin write.
        bin_data = rapidjson.loads(data)
        bin_data['start'] = int(bin_data['start'])
        bin_data['end'] = int(bin_data['end'])

        for line in ifile:
            chrom, start, stop, data = line.rstrip().split('\t')
            data = rapidjson.loads(data)
            data['start'] = int(data['start'])
            data['end'] = int(data['end'])
            if (abs(bin_data['mean'] - data['mean']) > args.fluctuation_limit) or (abs(bin_data['median'] - data['median']) > args.fluctuation_limit):
                write(ofile, bin_data)
                bin_data = data
            else:
                bin_data['end'] = data['end']
        write(ofile, bin_data)
    pysam.tabix_index(args.out_coverage_file, seq_col = 0, start_col = 1, end_col = 2, force = True)
