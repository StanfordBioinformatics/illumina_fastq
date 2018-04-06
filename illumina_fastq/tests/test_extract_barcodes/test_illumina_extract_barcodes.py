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
import gzip

import illumina_fastq.tests.test_extract_barcodes as testdir
import illumina_fastq.utils as fastq_utils

R1 = os.path.join(testdir.DATA_DIR,"INPUT/reads1.fastq")
R2 = os.path.join(testdir.DATA_DIR,"INPUT/reads2.fastq")
R1_OUTFILE = os.path.join(testdir.OUT_DIR,"forward_extract_3.fastq")
R2_OUTFILE = os.path.join(testdir.OUT_DIR,"reverse_extract_3.fastq")

if not os.path.exists(R1):
    raise Exception("Missing test input file {}".format(R1))
if not os.path.exists(R2):
    raise Exception("Missing test input file {}".format(R2))

if os.path.exists(R1_OUTFILE):
    os.remove(R1_OUTFILE)
if os.path.exists(R2_OUTFILE):
    os.remove(R2_OUTFILE)

outdir = testdir.OUT_DIR
class TestBarcodeExtract(unittest.TestCase):
    def setUp(self):
        """
        The input file R1 has five records, two of which have the barcode GTTACACT+GGATCTCG. The other three recs
        all contain unique barcodes.
        """

    def test_GTTACACT_GGATCTCG(self):
        """
        Barcode GTTACACT+GGATCTCG appears twice in the inputs.
        """
        outfile_prefix = os.path.join(outdir, "paired_end")

        r1_output_name = outfile_prefix + "_GTTACACT-GGATCTCG_R1.fastq.gz"
        if os.path.exists(r1_output_name):
            os.remove(r1_output_name)
        r2_output_name = outfile_prefix + "_GTTACACT-GGATCTCG_R2.fastq.gz"
        if os.path.exists(r2_output_name):
            os.remove(r2_output_name)

        cmd = "extract_barcodes.py --r1 {R1} --r2 {R2} --outfile-prefix {prefix} -b GTTACACT+GGATCTCG".format(
            R1=R1, R2=R2, prefix=outfile_prefix)
        subprocess.check_call(cmd, shell=True)

        r1_output = gzip.open(r1_output_name).read().strip()
        r2_output = gzip.open(r2_output_name).read().strip()
        self.assertEquals(
            r1_output,
            "@COOPER:74:HFTH3BBXX:3:1101:1052:1033 1:N:0:GTTACACT+GGATCTCG\nNGCCA\n+\n#AAFF\n@COOPER:74:HFTH3BBXX:3:1101:29833:1033 1:N:0:GTTACACT+GGATCTCG\nNATCC\n+\n#AAAA")
        self.assertEquals(
            r2_output,
            "@COOPER:74:HFTH3BBXX:3:1101:1052:1033 2:N:0:GTTACACT+GGATCTCG\nNGCCA\n+\n#AAFF\n@COOPER:74:HFTH3BBXX:3:1101:29833:1033 2:N:0:GTTACACT+GGATCTCG\nNTGTA\n+\n#AAAF")


if __name__ == '__main__':
    unittest.main()
