###
#Nathaniel Watson
#nathankw@stanford.edu
#Created Jan, 2017
###

import sys
import os
import datetime
import gzip
#from memory_profiler import profile

#@profile #used for memory_profiler

def parseIlluminaFastqAttLine(attLine):
	#Illumina FASTQ Att line format (as of CASAVA 1.8 at least):
	#  @<instrument-name>:<run ID>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read number>:<is filtered>:<control number>:<barcode sequence>
	uid = attLine.strip()
	header = uid.lstrip("@").split(":")
	dico = {}
	dico["instrument"] = header[0]
	dico["runId"] = header[1]
	dico["flowcellId"] = header[2]
	dico["lane"] = header[3] 
	dico["tile"] = header[4]
	dico["xpos"] = header[5]
	ypos,readNumber = header[6].split()
	dico["ypos"] = ypos 
	dico["readNumber"] = readNumber
	dico["isFiltered"] = header[7]
	dico["control"] = header[8]
	dico["barcode"] = header[9]
	return dico 

def get_pairedend_read_id(read_id):
	"""
	Function : Given either a forward read or reverse read identifier, returns the corresponding paired-end read identifier.
	Args     : read_id - str. forward read or reverse read identifier. This should be the entire title line of a FASTQ record,
	                     minus any trailing whitespace.
	Returns  : str. The pairend-end read identifier (title line). 
	Example  : Setting read_id to "@COOPER:74:HFTH3BBXX:3:1101:29894:1033 1:N:0:NATGAATC+NGATCTCG" will return 
	               @COOPER:74:HFTH3BBXX:3:1101:29894:1033 2:N:0:NATGAATC+NGATCTCG
	"""
	part1, part2 = read_id.strip().split()
	if part2.startswith("1"):
	  part2 = part2.replace("1","2",1)
	elif part2.startswith("2"):
	  part2 = part2.replace("2","1",1)
	else:
	  raise Exception("Unknown read number in {title}".format(title=read_id))
	return part1 + " " + part2


def isForwardRead(seqid):
  if seqid.split()[1].startswith("1"):
    return True
  return False

def yield_recs(fastqFile,log=sys.stdout,barcodes=[]):
	log.write("Parsing " + fastqFile + "\n")
	log.flush()
	fastqFileExt = os.path.splitext(fastqFile)[-1]
	if fastqFileExt == ".gz":
		fh = gzip.open(fastqFile)
	else:	
		fh = open(fastqFile)
	data = []
	lookup = {}
	count = 0
	lineCount = 0
	for line in fh:
		lineCount += 1
		count += 1
		data.append(line.strip())
		if count == 4:
			barcode = data[0].rsplit(":",1)[-1]
			if barcode in barcodes or not barcodes:
					yield data
			count = 0
			data = []
		if lineCount % 1000000 == 0:
			#every million lines
			log.write(str(datetime.datetime.now()) + ":  " + str(lineCount) + "\n")
			log.flush()
	fh.close()
	log.write("Finished parsing " + fastqFile + "\n")
	log.flush()
