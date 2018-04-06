# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
# nathan.watson86@gmail.com
###

# Created Jan, 2017

import datetime
import gzip
import os
import pdb
import sys
#from memory_profiler import profile

import illumina_fastq.utils as utils


class _FastqParseIter:
    def __init__(self, fastqparse_i):
        """
        Args:
            fastqparse_i: A FastqParse instance.
        """
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
    """
    Parses the records in an Illumina FASTQ file and stores all records or only those having
    specific barcodes. The sequence ID, sequence, and quality strings of each FASTQ record are
    stored in a list of lists of the form

        [ ["seqidA", "ACGT","#AAF"], ["seqidB", "GGAT"," #AAA"] ... ]

    This list of lists is stored as self.data. A lookup table (dict) is also stored as
    self.lookup. It is of the form

        { "seqidA": indexA, "seqidB": indexB, ... }

    where an index gives the position in the list of the record with the given sequence ID.
    The sequence ID is stored as the entire title line of a FASTQ record, minus any peripheral whitespace.

    Also supports indexing the returned instance object using the header line of a given sequence, i.e.
    if @GADGET:77:HFNLTBBXX:8:1101:30462:1279 1:N:0:NNAGCA is the read ID of a record that is present in a FASTQ file
    named reads.fq, then the following returns True:

        data = FastqParse("reads.fq")
        data["@GADGET:77:HFNLTBBXX:8:1101:30462:1279 1:N:0:NNAGCA"] #returns True


    Args:
        fastq: `str`. Path to the FASTQ file to be parsed. Accepts uncompressed or gzip
            compressed with a .gz extension.
        log: A file handle for logging. Defaults to STDOUT.
        extract_barcodes: `list` of one or more barcodes to extract from the FASTQ file.
            If the barcode is duel-indexed, separate with a '+', i.e. 'ATCGGT+GCAGCT',
            as this is how it is written in the FASTQ file.
        sample_size: `int`. Indicates the number of records from the start of the FASTQ file to
            parse. A Falsy value (the default) means that the entire FASTQ file will be parsed.
    """

    SEQID_KEY = "seqid"
    SEQ_KEY = "seq"
    QUAL_KEY = "qual"

    # SEQID_IDX, SEQ_IDX, and QUAL_IDX store the index position of the read ID, sequence string, and quality string, respectively,
    # of a given sublist in the list self.data.
    SEQID_IDX = 0
    SEQ_IDX = 1
    QUAL_IDX = 2

    def __init__(self, fastq, log=sys.stdout, extract_barcodes=[], sample_size=False):
        self.fastqFile = fastq
        self.barcodes = extract_barcodes
        self.sample_size = sample_size
        self.log = log
        self._parse()  # sets self.data
        # sets self.lookup.

    def __iter__(self):
        return _FastqParseIter(self)

    def __getitem__(self, seqid):
        return self.isRecordPresent(seqid)

    def __len__(self):
        return len(self.lookup)

    @classmethod
    def getPairedendReadId(cls, read_id):
        """
        Given either a forward read or reverse read identifier, returns the corresponding paired-end
        read identifier.

        Args:
            read_id: `str`. The forward read or reverse read identifier. This should be the
                entire title line of a FASTQ record, minus any trailing whitespace.
        Returns:
            `str`: The paired-end read identifier (title line).
        Example:
            Setting read_id to "@COOPER:74:HFTH3BBXX:3:1101:29894:1033 1:N:0:NATGAATC+NGATCTCG" will return
            @COOPER:74:HFTH3BBXX:3:1101:29894:1033 2:N:0:NATGAATC+NGATCTCG
        """
        return utils.getPairedendReadId(read_id)

    @classmethod
    def formatRecordForOutput(cls, record):
        return "\n".join([record[FastqParse.SEQID_KEY], record[FastqParse.SEQ_KEY],
                          "+", record[FastqParse.QUAL_KEY]]) + "\n"

    def printRecord(self, seqid, outfh):
        outfh.write(FastqParse.formatRecordForOutput(self.getRecord(seqid)))

    @classmethod
    def parseIlluminaFastqAttLine(cls, attLine):
        # Illumina FASTQ Att line format (as of CASAVA 1.8 at least):
        #  @<instrument-name>:<run ID>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read number>:<is filtered>:<control number>:<barcode sequence>
        return utils.parseIlluminaFastqAttLine(attLine)

    # @profile #used for memory_profiler
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
                barcode = line.rsplit(":", 1)[-1]
                #self.data[uid] = {"name": line}
                #self.data[uid] = {}
            elif count == 2:
                seq = line
                #self.data[uid]["seq"] = line
            elif count == 4:
                #self.data[uid]["qual"] = line
                if barcode in self.barcodes or not self.barcodes:
                    self.data.append([uid, seq, line])
                    hash_id = hash(uid)
                    if hash_id in self.lookup:
                        raise Exception(
                            "Found multiple entires for the FASTQ record having title line '{uid}'.".format(
                                uid=uid))
                    self.lookup[hash_id] = len(self.data) - 1
                    if self.sample_size and len(self.lookup) == self.sample_size:
                        break
                count = 0
            if lineCount % 1000000 == 0:
                # every million lines
                self.log.write(str(datetime.datetime.now()) + ":  " + str(lineCount) + "\n")
                self.log.flush()
        fh.close()
        self.log.write("Finished parsing " + self.fastqFile + "\n")
        self.log.flush()

    def isRecordPresent(self, title_line):
        if hash(title_line) in self.lookup:
            return True
        return False

    def getRecord(self, title_line):
        rec = self.data[self.lookup[hash(title_line)]]
        return self._formatRecord(rec)

    @classmethod
    def isForwardRead(cls, seqid):
        return utils.isForwardRead(seqid)

    def _formatRecord(self, rec):
        """
        Args : rec - A sublist from self.data.
        """
        return {FastqParse.SEQID_KEY: rec[FastqParse.SEQID_IDX],
                FastqParse.SEQ_KEY: rec[FastqParse.SEQ_IDX], FastqParse.QUAL_KEY: rec[FastqParse.QUAL_IDX]}

    def barcodeHist(self):
        bcDico = {}
        count = 0
        for rec in self:
            att_line = FastqParse.parseIlluminaFastqAttLine(rec[FastqParse.SEQID_KEY])
            barcode = att_line["barcode"]
            if barcode not in bcDico:
                bcDico[barcode] = 0
            bcDico[barcode] += 1
        return bcDico
