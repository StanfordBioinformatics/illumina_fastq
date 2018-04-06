# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
# nathan.watson86@gmail.com
###

# Created Jan, 2017

import sys
import os
import datetime
import gzip
#from memory_profiler import profile

# @profile #used for memory_profiler


def parseIlluminaFastqAttLine(attLine):
    """
    Given the title line of a FASTQ record, tonizes the line and stores the tokens in a dict.
    The Illumina FASTQ Att line format (as of CASAVA 1.8 at least) is:

      @<instrument-name>:<run ID>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read number>:<is filtered>:<control number>:<barcode sequence>

    Args:
        attLine: `str`. The title line of a FASTQ record, minus any trailing whitespace.
    Returns:
        `dict`. The keys are:

        1. instrument,
        2. runId,
        3. flowcellId,
        4. lane,
        5. tile,
        6. xpos,
        7. ypos,
        8. readNumber,
        9. isFiltered,
       10. control,
       11. barcode

    """
    uid = attLine.strip()
    header = uid.lstrip("@").split(":")
    dico = {}
    dico["instrument"] = header[0]
    dico["runId"] = header[1]
    dico["flowcellId"] = header[2]
    dico["lane"] = header[3]
    dico["tile"] = header[4]
    dico["xpos"] = header[5]
    ypos, readNumber = header[6].split()
    dico["ypos"] = ypos
    dico["readNumber"] = readNumber
    dico["isFiltered"] = header[7]
    dico["control"] = header[8]
    dico["barcode"] = header[9]
    return dico


def getPairedendReadId(read_id):
    """
    Given either a forward read or reverse read identifier, returns the corresponding paired-end
    read identifier.

    Args:
        read_id: `str`. The forward read or reverse read identifier. This should be the entire
            title line of a FASTQ record, minus any trailing whitespace.
    Returns:
        `str`: The pairend-end read identifier (title line).

    Example:
        Setting read_id to "@COOPER:74:HFTH3BBXX:3:1101:29894:1033 1:N:0:NATGAATC+NGATCTCG" will return
        @COOPER:74:HFTH3BBXX:3:1101:29894:1033 2:N:0:NATGAATC+NGATCTCG
    """
    part1, part2 = read_id.strip().split()
    if part2.startswith("1"):
        part2 = part2.replace("1", "2", 1)
    elif part2.startswith("2"):
        part2 = part2.replace("2", "1", 1)
    else:
        raise Exception("Unknown read number in {title}".format(title=read_id))
    return part1 + " " + part2


def isForwardRead(seqid):
    """Indicates whether the passed-in read identifier is a forward or reverse read identifier.

    Args:
        seqid: `str`. A read identifier of a FASTQ record.

    Returns:
        `bool`: True if a forward read identifier, False otherwise.
    """
    if seqid.split()[1].startswith("1"):
        return True
    return False


def yieldRecs(fastqFile, log=sys.stdout, barcodes=[]):
    """
    A generator function that reads a FASTQ file and yields records, one at a time.
    The records to yield can be restricted to the specified set of barcodes.

    Args:
        fastqFile: `str`. Path to the FASTQ file to parse.
        log: A file handle to write log messages to. Defaults to STDOUT.

    Yields:
        A list containing one element per line of a FASTQ record. Each element is whitespace stripped.
    """
    log.write("Parsing " + fastqFile + "\n")
    log.flush()
    fastqFileExt = os.path.splitext(fastqFile)[-1]
    if fastqFileExt == ".gz":
        fh = gzip.open(fastqFile)
    else:
        fh = open(fastqFile)
    data = []
    count = 0
    lineCount = 0
    for line in fh:
        lineCount += 1
        count += 1
        data.append(line.strip())
        if count == 4:
            barcode = data[0].rsplit(":", 1)[-1]
            if barcode in barcodes or not barcodes:
                yield data
            count = 0
            data = []
        if lineCount % 1000000 == 0:
            # every million lines
            log.write(str(datetime.datetime.now()) + ":  " + str(lineCount) + "\n")
            log.flush()
    fh.close()
    log.write("Finished parsing " + fastqFile + "\n")
    log.flush()
