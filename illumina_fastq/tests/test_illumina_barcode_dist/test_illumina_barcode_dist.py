#!/bin/env python
# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
# nathan.watson86@gmail.com
###

# Created Oct 9, 2017

import subprocess
import os
import unittest

import illumina_fastq.tests.test_illumina_barcode_dist as testdir
from illumina_fastq.illumina_fastq_parse import FastqParse

R1 = "INPUT/reads1.fastq"
OUTDIR = testdir.OUT_DIR
OUTFILE = os.path.join(OUTDIR,"bcdist.txt")

if not os.path.exists(R1):
    raise Exception("Missing test input file {}".format(R1))


class TestBarcodeDist(unittest.TestCase):
    def setUp(self):
        """
        The input file R1 has five records, two of which have the barcode GTTACACT+GGATCTCG. The other three recs
        all contain unique barcodes.
        """
        if os.path.exists(OUTFILE):
            os.remove(OUTFILE)

        self.hist = FastqParse(fastq=R1).barcodeHist()
        #cmd = "illuminaBarcodeDist.py -i {infile} -o {outfile} ".format(infile=R1,outfile=OUTFILE)
        # subprocess.check_call(cmd,shell=True)

    def test_correct_dist_(self):
        res = {
            "GTTACACT+GGATCTCG": 2,
            "NATGAATC+NGATCTCG": 1,
            "CATGAATC+TGATCTCG": 1,
            "CATGAATC+GGATCTCG": 1}
        self.assertEqual(self.hist, res)


if __name__ == '__main__':
    unittest.main()
