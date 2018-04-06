#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
# nathan.watson86@gmail.com
###

# Created Feb 7, 2014

import sys
import operator
import argparse
import os
import glob
import subprocess
import gzip

from illumina_fastq.illumina_fastq_parse import FastqParse

DEFAULT_SAMPLE_SIZE = 100000

"""
Given a FASTQ file containing unmatched reads, tabulates the frequencies at which each barcode is 
present. The first line of each FASTQ record must be in the standard Illumina format:
	@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>[:<UMI>] <read>:<is filtered>:<control number>:<index sequence>
where anything in [] is optional. The output file will have 3 tab-delimite fields being
	1) Barcode
	2) Freq
	3) Relative_Freq%

|
"""

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-i',
        '--infile',
        required=True,
        help="Input FASTQ file. May be gzip'd with a .gz extension.")
    parser.add_argument(
        '-o',
        '--outfile',
        help="Output file name with barcode histogram. Default is --infile name with the added extension '_barcodeHist.txt'.")
    parser.add_argument(
        "-s",
        "--sample-size",
        type=int,
        default=DEFAULT_SAMPLE_SIZE,
        help="The number of reads to use to create the distribution, taken from the start of the file. Default is {:,}".format(DEFAULT_SAMPLE_SIZE))

    parser.add_argument('-v', '--version', action="version", version='%(prog)s 0.2')

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    
    fqfile = args.infile
    outfile = args.outfile
    if not outfile:
        outfile = fqfile + "_barcodeHist.txt"
    sample_size = args.sample_size
    
    hist = FastqParse(fastq=fqfile, sample_size=sample_size).barcodeHist()
    num_recs = sum(hist.values())
    
    fout = open(outfile, 'w')
    outputHeader = ["Barcode", "Freq", "Relative_Freq%"]
    fout.write("\nTotal number of records sampled: {0}\n".format(num_recs))
    fout.write("\n")
    fout.write("\t".join(outputHeader) + "\n")
    for bc, cnt in sorted(hist.items(), key=operator.itemgetter(1), reverse=True):
        perc = cnt / float(num_recs)
        fout.write("{bc}\t{cnt}\t{perc:.2%}\n".format(bc=bc, cnt=cnt, perc=perc))
    fout.close()

if __name__ == "__main__":
    main()
