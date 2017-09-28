#!/usr/bin/env python

###
#Nathaniel Watson
#nathankw@stanford.edu
#Stanford Center for Genetics and Personalized Medicine
#2017-09-27
###

import os
import sys
import datetime
import gzip
from argparse import ArgumentParser

from illumina_fastq.illumina_fastq_parse import FastqParse

FASTQ_EXT =  ".fastq"
GZIP_EXT = ".gz"
R1 = "R1"
R2 = "R2"

description ="Say you have a foward read and you want to fetch its corresponding paired-end read. This script does just that, in batch. You provide a FASTQ file containing a subset of forward and reverse reads (or either thereof), and this script extracts the corresponding mates from the parent FASTQ files. A new pair of FASTQ files will be written (forwards reads and reverse reads files), sorted by query-name. This only works for Illumina formatted reads (with regard to the header lines)."

parser = ArgumentParser(description=description)
parser.add_argument("-q","--query-reads",required=True,help="The FASTQ file containing forward, reverse, or a mix of forward and reverse reads whose mates need to be also fetched.")
parser.add_argument("-f","--forward-reads",required=True,help="The parent foward reads FASTQ file")
parser.add_argument("-r","--reverse-reads",required=True,help="The parent reverse reads FASTQ file")
parser.add_argument("--outfile-f",required=True,help="Name of the output file for the foward extracted READS (includes any forward reads present in the query set).")
parser.add_argument("--outfile-r",required=True,help="Name of the output file for the reverse extracted READS (includes any reverse reads present in the query set).")
args = parser.parse_args()

q_file = args.query_reads
r1_file = args.forward_reads
r2_file = args.reverse_reads

query = FastqParse(q_file)
r1 = FastqParse(r1_file)
r2 = FastqParse(r2_file)


fout_f = open(args.outfile_f,'w')
fout_r = open(args.outfile_r,'w')

for rec in query:
	seqid = rec[FastqParse.SEQID_KEY]
	mateid = query.get_pairedend_read_id(seqid)
  #check if mate is already in the query file
	if query[mateid]:
		continue
	if query.isForwardRead(seqid):
		query.printRecord(seqid,fout_f)
		r2.printRecord(mateid,fout_r)
	else:
		query.printRecord(seqid,fout_r)
		r1.printRecord(mateid,fout_f)

fout_f.close()
fout_r.close()
