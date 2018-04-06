#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
# nathan.watson86@gmail.com
###

# Created Sep, 2017

import os
import sys
import datetime
import gzip
import subprocess
from argparse import ArgumentParser

from illumina_fastq.illumina_fastq_parse import FastqParse

description = "Generates a file containing the FASTQ record identifiers (one per line) from the passed in FASTQ files."

parser = ArgumentParser(description=description)
parser.add_argument(
    "-i",
    "--infiles",
    nargs="+",
    required=True,
    help="One or more FASTQ input files. Can be gzip-compressed with a .gz extension.")
parser.add_argument(
    "-o",
    "--outfile",
    required=True,
    help="The output file containing a single column, where each row is a FASTQ identifier.")
args = parser.parse_args()

infiles = args.infiles
outfile = args.outfile
fout = open(outfile, 'w')

for fqfile in args.infiles:
    records = FastqParse(fastq=fqfile)
    for rec in records:
        fout.write(rec[FastqParse.SEQID_KEY] + "\n")
fout.close()
