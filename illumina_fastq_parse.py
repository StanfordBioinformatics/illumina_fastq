
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



class _FastqParseIter:
	def __init__(self,fastqparse_i):
		self.idx = -1
		self.fpi = fastqparse_i
		self.seqids = self.fpi.lookup.keys()
		self.len_seqids = len(self.seqids)

	def next(self):	
		self.idx += 1
		if self.idx == self.len_seqids:
			raise StopIteration
		seqid_hashval = self.seqids[self.idx]
		seqid_idx = self.fpi.lookup[seqid_hashval]
		return self.fpi._formatRecord(self.fpi.data[seqid_idx])


class FastqParse:

	SEQID_KEY = "seqid"
	SEQ_KEY = "seq"
	QUAL_KEY = "qual"

	#SEQID_IDX, SEQ_IDX, and QUAL_IDX store the index position of the read ID, sequence string, and quality string, respectively,
	# of a given sublist in the list self.data.
	SEQID_IDX = 0 
	SEQ_IDX = 1
	QUAL_IDX = 2

	def __init__(self,fastq,log=sys.stdout,extract_barcodes=[]):
		"""
		Function : Parses the records in an Illumina FASTQ file and stores all records or only those having specific barcodes.
							 The sequence ID, sequence, and quality strings of each FASTQ record are stored in a list of lists of the form
																		[ ["seqidA", "ACGT","#AAF"], ["seqidB", "GGAT"," #AAA"] ... ] 
							 This list of lists is stored as self.data. A lookup table (dict) is also stored as self.lookup. It is of the form 
																		{ "seqid_x": index_x, "seqid_y": index_y, ... }
							 where an index gives the position in the list of the record with the given sequence ID. The sequence ID is stored
							 as the entire title line of a FASTQ record.

               Also supports indexing the returned instance object using the header line of a given sequence, i.e. 
               if @GADGET:77:HFNLTBBXX:8:1101:30462:1279 1:N:0:NNAGCA is the read ID of a record that is present in a FASTQ file 
               named reads.fq, then the following returns True:
                   
                   data = FastqParse("reads.fq")
                   data["@GADGET:77:HFNLTBBXX:8:1101:30462:1279 1:N:0:NNAGCA"] #returns True
                   

		Args     : fastq - The FASTQ file to be parsed. Accepts uncompressed or gzip compressed with a .gz extension.
							 log - file handle for logging. Defaults to sys.stdout.
							 extract_barcodes - list of one or more barcodes to extract from the FASTQ file. If the barcode is duel-indexed, separate
							     them with a '+', i.e. 'ATCGGT+GCAGCT', as this is how it is written in the FASTQ file. 
		"""
		self.fastqFile = fastq
		self.barcodes = extract_barcodes
		self.log = log
		self._parse() #sets self.data
									#sets self.lookup.

	def __iter__(self):
		return _FastqParseIter(self)

	@classmethod
	def get_pairedend_read_id(cls,read_id):
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

	@classmethod
	def formatRecordForOutput(cls,record):
		return "\n".join([record[FastqParse.SEQID_KEY],record[FastqParse.SEQ_KEY],"+",record[FastqParse.QUAL_KEY]]) + "\n"

	def printRecord(self,seqid,outfh):
		outfh.write(FastqParse.formatRecordForOutput(self.getRecord(seqid)))

	@classmethod 
	def parseIlluminaFastqAttLine(cls,attLine):
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

	#@profile #used for memory_profiler
	def _parse(self):
		self.log.write("Parsing " + self.fastqFile + "\n")
		self.log.flush()
		fastqFileExt = os.path.splitext(self.fastqFile)[-1]
		if fastqFileExt == ".gz":
			fh = gzip.open(self.fastqFile)
		else:	
			fh = open(self.fastqFile)
		self.data = []
		self.lookup = {}
		count = 0
		lineCount = 0
		for line in fh:
			lineCount += 1
			count += 1
			line = line.strip()
			if count == 1:
				#uid = lineCount
				uid = line
				barcode = line.rsplit(":",1)[-1]
				#self.data[uid] = {"name": line}
				#self.data[uid] = {}
			elif count == 2:
				seq = line
				#self.data[uid]["seq"] = line
			elif count == 4:
				#self.data[uid]["qual"] = line
				if barcode in self.barcodes or not self.barcodes:
						self.data.append([uid,seq,line])
						hash_id = hash(uid)
						if hash_id in self.lookup:
							raise Exception("Found multiple entires for the FASTQ record having title line '{uid}'.".format(uid=uid))
						self.lookup[hash_id] = len(self.data) - 1
				count = 0
			if lineCount % 1000000 == 0:
				self.log.write(str(datetime.datetime.now()) + ":  " + str(lineCount) + "\n")
				self.log.flush()
		fh.close()
		self.log.write("Finished parsing " + self.fastqFile + "\n")
		self.log.flush()
	
	def __getitem__(self,seqid):
		try:
			self.getRecord(seqid)
		except KeyError:
			return False
		return True
		

	def getRecord(self,title_line):
		rec = self.data[self.lookup[hash(title_line)]]
		return self._formatRecord(rec)

	@classmethod
	def isForwardRead(cls,seqid):
		if seqid.split()[1].startswith("1"):
			return True
		return False

	def _formatRecord(self,rec):
		"""
		Args : rec - A sublist from self.data.
		"""
		return { FastqParse.SEQID_KEY: rec[FastqParse.SEQID_IDX], FastqParse.SEQ_KEY: rec[FastqParse.SEQ_IDX], FastqParse.QUAL_KEY: rec[FastqParse.QUAL_IDX] }
	
#Total number of lines in SCGPM_MD-DNA-1_HFTH3_L3_unmatched_R1.fastq is 347,060,820.


if __name__ == "__main__":
  FastqParse(sys.argv[1])
