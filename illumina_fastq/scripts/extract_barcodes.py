#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
# nathan.watson86@gmail.com
###

# Created Jan, 2016

"""
Extracts FASTQ records matching the specified barcodes from a pair of unmatched read FASTQ files. 
For each specified barcode, records matching the barcode will be written to a new pair of FASTQ 
files. All output files are compressed with gzip. It is expected that both input FASTQ files 
(the forward and reverse read files) are sorted by read name and that their aren't any missing 
mates, otherwise an Exception will be raised.

|
"""

import os
import sys
import datetime
import gzip
import argparse

import illumina_fastq.utils as fastq_utils

GZIP_FASTQ_EXT = ".fastq.gz"
R1 = "R1"
R2 = "R2"

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--r1", required=True, help="""
    FASTQ file containing the (forward) reads.""")

    parser.add_argument("--r2", required=True, help="""
    FASTQ file containing the reverse reads.""")

    parser.add_argument("--outdir", help="""
    The pre-existing directory to output the FASTQ files containing the extracted barcodes. Defaults to the current working directory.""")

    parser.add_argument("--outfile-prefix", required=True, help="""
    The file prefix of each output FASTQ file for a given barcode. The barcode name will be 
    appended to this prefix, as well as the read number. For example, setting the outfile prefix 
    to 'output' would result in 'output_${barcode}_R1.fastq' and 'output_${barcode}_R2.fastq'.""")

    parser.add_argument("-b", "--barcodes", nargs="+", help="""
    One or more barcodes to extract from the input FASTQ file(s).""")

    return parser

def main():
    parser = get_parser()                                                                              
    args = parser.parse_args()

    barcodes = args.barcodes
    r1_file = args.r1
    r2_file = args.r2
    outdir = args.outdir
    if not outdir:
        outdir = os.getcwd()
    else:
        if not os.path.exists(outdir):
            raise Exception("The path {path} passed to --outdir does not exist!".format(path=outdir))
    
    outfile_prefix = args.outfile_prefix
    barcodes = [x.replace("-", "+") for x in barcodes]
    # In the title line of Illumina FASTQ records, a duel-indexed barcode is
    # separated with a '+', i.e. ATC+CGA.
    
    if not os.path.exists(r1_file):
        raise Exception("{r1_file} provided to --r1 doesn't exist!".format(r1_file=r1_file))
    if not os.path.exists(r2_file):
        raise Exception("{r2_file} proided to --r2 doens't exist!".format(r2_file=r2_file))
    
    start_time = datetime.datetime.now()
    
    r1_records = fastq_utils.yieldRecs(fastqFile=r1_file, barcodes=barcodes)
    r2_records = fastq_utils.yieldRecs(fastqFile=r2_file, barcodes=barcodes)
    
    file_handles = {}
    for barcode in barcodes:
        file_handles[barcode] = {}
        outfile_name = os.path.join(outdir, outfile_prefix + "_" + barcode.replace("+", "-"))
        outfile_name += "_" + R1 + GZIP_FASTQ_EXT
        file_handles[barcode][R1] = gzip.open(outfile_name, "wb")
        file_handles[barcode][R2] = gzip.open(outfile_name.replace(R1, R2), "wb")
    
    output_barcode_counts = {}
    for barcode in barcodes:
        output_barcode_counts[barcode] = 0
    
    print("Extracting barcodes of interest")
    sys.stdout.flush()
    
    while True:
        try:
            f_rec = r1_records.next()  # forward read record
            r_rec = r2_records.next()  # reverse read record
        except StopIteration:
            break
        f_id = f_rec[0]
        if fastq_utils.getPairedendReadId(f_id) != r_rec[0]:
            raise Exception(
                "Input FASTQ files out of order. They need to be sorted by name and all pairs should be present.")
        barcode = fastq_utils.parseIlluminaFastqAttLine(f_id)["barcode"]
    
        file_handles[barcode][R1].write("\n".join(f_rec) + "\n")
        file_handles[barcode][R2].write("\n".join(r_rec) + "\n")
        output_barcode_counts[barcode] += 1
    
    to_compress = []
    for barcode in file_handles:
        for output_key in file_handles[barcode]:
            handle = file_handles[barcode][output_key]
            handle.close()
            if os.stat(handle.name).st_size == 0:
                os.remove(handle.name)
            else:
                to_compress.append(handle.name)
    
    end_time = datetime.datetime.now()
    print("\n")
    sys.stdout.flush()
    print("Elapsed time: {time}\n".format(time=str(end_time - start_time)))
    sys.stdout.flush()
    print("Output Statistics")
    sys.stdout.flush()
    for barcode in output_barcode_counts:
        print(barcode + ": " + str(output_barcode_counts[barcode]))
    print("\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
