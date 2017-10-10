#!/usr/bin/env python

###
#Nathaniel Watson
#nathankw@stanford.edu
#Stanford Center for Genetics and Personalized Medicine
#2016-01-06
###

import os
import sys
import datetime
import gzip
import subprocess
from argparse import ArgumentParser

import illumina_fastq.utils as fastq_utils

description ="Extracts FASTQ records matching the specified barcodes from a pair of unmatched read FASTQ files. For each specified barcode, records matching the barcode will be written to a new pair of FASTQ files. All output files are compressed with gzip. If one of the reads in a pair that matches a given barcode isn't present in the unmatched reads, then an exception is raised. It is expected that both input FASTQ files (the forward and reverse read files) are sorted by read name and that their aren't any missing mates."

parser = ArgumentParser(description=description)
parser.add_argument("--r1",required=True,help="FASTQ file containing the (forward) reads.")
parser.add_argument("--r2",required=True,help="FASTQ file containing the reverse reads.")
parser.add_argument("--outdir",help="The pre-existing directory to output the FASTQ files containing the extracted barcodes. Defaults to the current working directory.")
parser.add_argument("--outfile-prefix",required=True,help="The file prefix of each output FASTQ file for a given barcode. The barcode name will be appended to this prefix, as well as the read number. For example, setting the outfile prefix to 'output' would result in 'output_${barcode}_R1.fastq' and 'output_${barcode}_R2.fastq'.")
parser.add_argument("-b","--barcodes",nargs="+",help="One or more barcodes to extract from the input FASTQ file(s).")

FASTQ_EXT =  ".fastq"
R1 = "R1"
R2 = "R2"

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
barcodes = [x.replace("-","+") for x in barcodes] 
#In the title line of Illumina FASTQ records, a duel-indexed barcode is separated with a '+', i.e. ATC+CGA.

if not os.path.exists(r1_file):
	raise Exception("{r1_file} provided to --r1 doesn't exist!".format(r1_file=r1_file))
if not os.path.exists(r2_file):
	raise Exception("{r2_file} proided to --r2 doens't exist!".format(r2_file=r2_file))

start_time = datetime.datetime.now()

r1_records = fastq_utils.yield_recs(fastqFile=r1_file,barcodes=barcodes)
r2_records = fastq_utils.yield_recs(fastqFile=r2_file,barcodes=barcodes)

file_handles = {}
for barcode in barcodes:
	file_handles[barcode] = {}
	outfile_name = os.path.join(outdir,outfile_prefix + "_" + barcode.replace("+","-"))
	outfile_name += "_" + R1 + FASTQ_EXT
	file_handles[barcode][R1] = open(outfile_name,"w")
	file_handles[barcode][R2] = open(outfile_name.replace(R1,R2),"w")

output_barcode_counts = {}
for barcode in barcodes:
	output_barcode_counts[barcode] = 0

print("Extracting barcodes of interest")
sys.stdout.flush()

while True:
	try:
		f_rec = r1_records.next() #forward read record
		r_rec = r2_records.next() #reverse read record
	except StopIteration:
		break
	f_id = f_rec[0]
	if fastq_utils.get_pairedend_read_id(f_id) != r_rec[0]:
		raise Exception("Input FASTQ files out of order. They need to be sorted by name and all pairs should be present.")
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

print("Compressing output files with gzip")
sys.stdout.flush()
for i in to_compress:
	cmd = "gzip {i}".format(i=i)
	print("Compressing {i} with command '{cmd}'.".format(i=i,cmd=cmd))
	sys.stdout.flush()
	popen = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	stdout,stderr = popen.communicate()
	retcode = popen.returncode
	if retcode:
		raise Exception("Failed to compress {i}. Stderr is '{stderr}'.".format(i=i,stderr=stderr))	
