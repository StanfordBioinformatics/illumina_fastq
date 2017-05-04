#!/usr/bin/env python

###AUTHOR###
#Nathaniel Watson
###DATE###
#Feb. 7, 2014


import sys
import operator
import argparse
import os
import glob
import subprocess
import gzip
#from gbsc import confParse 

#config = "/srv/gs1/software/conf/gbsc/global.txt"
#conf = confParse.Parse(config)
#pubDir = conf.getValue("illuminaPublished")

def barcodeHist(fqFile,outfile,sample_size):
	if fqfile.endswith(".gz"):
		fh = gzip.open(fqFile,'r')
	else:
		fh = open(fqFile,'r')
	fout = open(outfile,'w')
	outputHeader = ["Barcode","Freq","Relative_Freq%"]
	bcDico = {}
	count = 0
	record_count = 0
	for line in fh:
		if record_count >= sample_size:
			break
		line = str(line) #is a bytes object if opened with gzip in Python3
		line = line.strip()
		if not line:
			continue
		count += 1
		if count == 4:
			record_count += 1
			count = 0
			continue
		if count in [2,3]:
			continue
		barcode = line.split(":")[-1]
		if barcode not in bcDico:
			bcDico[barcode] = 0
		bcDico[barcode] += 1
	fh.close()
	numRecs = float(sum(bcDico.values()))

	fout.write("\nTotal number of records sampled: {0}\n".format(int(numRecs)))
	fout.write("\n")
	fout.write("\t".join(outputHeader) + "\n")
	for index,cnt in sorted(bcDico.items(),key=operator.itemgetter(1),reverse=True):
		perc = cnt/numRecs
		fout.write("{index}\t{cnt}\t{perc:.2%}\n".format(index=index,cnt=cnt,perc=perc))
	fout.close()

description="Given the full path to a run name, tabulates the frequencies at which each barcode in the Illumina Unmatched FASTQ files (from demultiplexing) are present. The first line of a FASTQ record) must be formatted as @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>."
parser = argparse.ArgumentParser(description=description)
parser.add_argument('-i','--infile',required=True,help="Input FASTQ file.")
parser.add_argument('-o','--outfile',help="Output file with barcode histogram. Default is --infile name with the added extension '_barcodeHist.txt'.")
parser.add_argument("-s","--sample-size",type=int,default=100000,help="The number of reads to use to create the distribution, taken from the start of the file. Default is %(default)s.")
parser.add_argument('-v','--version',action="version",version='%(prog)s 0.2')
args = parser.parse_args()

fqfile = args.infile
outfile = args.outfile
if not outfile:
	outfile = fqfile + "_barcodeHist.txt"
sample_size = args.sample_size

barcodeHist(fqfile,outfile,sample_size)

