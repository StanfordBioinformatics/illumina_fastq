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

from illumina_fastq.illumina_fastq_parse import FastqParse

R1 = "INPUT/reads1.fastq"
R2 = "INPUT/reads2.fastq"

if not os.path.exists(R1):
    raise Exception("Missing test input file {}".format(R1))
if not os.path.exists(R2):
    raise Exception("Missing test input file {}".format(R2))


class TestOutputContainsQueries(unittest.TestCase):
    def setUp(self):
        """
        The query file has the two records below from the reads1.fastq file:
                1) @COOPER:74:HFTH3BBXX:3:1101:1052:1033 1:N:0:NTTACACT+NGATCTCG
                2) @COOPER:74:HFTH3BBXX:3:1101:29833:1033 1:N:0:NAGATTCC+NAGCTATA
        Both R1 and R2 FASTQ files contain the respective paired-end reads entries. The purpose
        of this test class is to ensure that both records are extracted to the output files.
        """
        query = "query_forward.fastq"
        self.query_fid1 = "@COOPER:74:HFTH3BBXX:3:1101:1052:1033 1:N:0:NTTACACT+NGATCTCG"
        self.query_fid2 = "@COOPER:74:HFTH3BBXX:3:1101:29833:1033 1:N:0:NAGATTCC+NAGCTATA"
        self.query_rid1 = "@COOPER:74:HFTH3BBXX:3:1101:1052:1033 2:N:0:NTTACACT+NGATCTCG"
        self.query_rid2 = "@COOPER:74:HFTH3BBXX:3:1101:29833:1033 2:N:0:NAGATTCC+NAGCTATA"
        self.fout_name = "forward_extract.fastq"
        self.rout_name = "reverse_extract.fastq"
        if os.path.exists(self.fout_name):
            os.remove(self.fout_name)
        if os.path.exists(self.rout_name):
            os.remove(self.rout_name)
        cmd = "add_mate_reads.py -q {q} -f {f} -r {r} --outfile-f {fout_name} --outfile-r {rout_name}".format(
            q=query, f=R1, r=R2, fout_name=self.fout_name, rout_name=self.rout_name)
        print("Running command: {}").format(cmd)
        subprocess.check_call(cmd, shell=True)

    def test_outfile_identifiers(self):
        r1 = FastqParse(self.fout_name)
        r2 = FastqParse(self.rout_name)
        self.assertEqual(r1[self.query_fid1], True)
        self.assertEqual(r1[self.query_fid2], True)
        self.assertEqual(r2[self.query_rid1], True)
        self.assertEqual(r2[self.query_rid2], True)

    def test_output_num_recs(self):
        r1 = FastqParse(self.fout_name)
        r2 = FastqParse(self.rout_name)
        self.assertEqual(len(r1), 2, "Incorrect number of records in forward reads output file.")
        self.assertEqual(len(r2), 2, "Incorrect number of records in reverse reads output file.")


class TestOutputContainsQueries_2(unittest.TestCase):
    def setUp(self):
        """
        Same as the class above (TestOutputContainsQueries) except that this time we're testing
        that we get the same results when starting from two reverse reads in the query instead of the
        two forward reads.
        The query file has the two records below from the reads2.fastq file:
                1) @COOPER:74:HFTH3BBXX:3:1101:1052:1033 2:N:0:NTTACACT+NGATCTCG
                2) @COOPER:74:HFTH3BBXX:3:1101:29833:1033 2:N:0:NAGATTCC+NAGCTATA
        Both R1 and R2 FASTQ files containt the respective paired-end reads entries. The purpose
        of this test class is to ensure that both records are extracted to the output files.
        """
        query = "query_reverse.fastq"
        self.query_fid1 = "@COOPER:74:HFTH3BBXX:3:1101:1052:1033 1:N:0:NTTACACT+NGATCTCG"
        self.query_fid2 = "@COOPER:74:HFTH3BBXX:3:1101:29833:1033 1:N:0:NAGATTCC+NAGCTATA"
        self.query_rid1 = "@COOPER:74:HFTH3BBXX:3:1101:1052:1033 2:N:0:NTTACACT+NGATCTCG"
        self.query_rid2 = "@COOPER:74:HFTH3BBXX:3:1101:29833:1033 2:N:0:NAGATTCC+NAGCTATA"
        self.fout_name = "forward_extract_2.fastq"
        self.rout_name = "reverse_extract_2.fastq"
        if os.path.exists(self.fout_name):
            os.remove(self.fout_name)
        if os.path.exists(self.rout_name):
            os.remove(self.rout_name)
        cmd = "add_mate_reads.py -q {q} -f {f} -r {r} --outfile-f {fout_name} --outfile-r {rout_name}".format(
            q=query, f=R1, r=R2, fout_name=self.fout_name, rout_name=self.rout_name)
        print("Running command: {}").format(cmd)
        subprocess.check_call(cmd, shell=True)

    def test_outfile_identifiers(self):
        r1 = FastqParse(self.fout_name)
        r2 = FastqParse(self.rout_name)
        self.assertEqual(r1[self.query_fid1], True)
        self.assertEqual(r1[self.query_fid2], True)
        self.assertEqual(r2[self.query_rid1], True)
        self.assertEqual(r2[self.query_rid2], True)

    def test_output_num_recs(self):
        r1 = FastqParse(self.fout_name)
        r2 = FastqParse(self.rout_name)
        self.assertEqual(len(r1), 2, "Incorrect number of records in forward reads output file.")
        self.assertEqual(len(r2), 2, "Incorrect number of records in reverse reads output file.")


class TestOutputContainsQueries_3(unittest.TestCase):
    def setUp(self):
        """
        Here, one of the queries is present in the parent FASTQ files, and one of them isn't.
        The query file has the record below from the reads2.fastq file:
                1) @COOPER:74:HFTH3BBXX:3:1101:29833:1033 2:N:0:NAGATTCC+NAGCTATA
        And the non-existant record:
                2) @not_present_
        Need to make sure that the output files contain only the record from the first read.
        """
        query = "query_notpresent.fastq"
        self.query_fid1 = "@COOPER:74:HFTH3BBXX:3:1101:29833:1033 1:N:0:NAGATTCC+NAGCTATA"
        self.query_rid1 = "@COOPER:74:HFTH3BBXX:3:1101:29833:1033 2:N:0:NAGATTCC+NAGCTATA"
        self.fout_name = "forward_extract_3.fastq"
        self.rout_name = "reverse_extract_3.fastq"
        if os.path.exists(self.fout_name):
            os.remove(self.fout_name)
        if os.path.exists(self.rout_name):
            os.remove(self.rout_name)
        cmd = "add_mate_reads.py -q {q} -f {f} -r {r} --outfile-f {fout_name} --outfile-r {rout_name}".format(
            q=query, f=R1, r=R2, fout_name=self.fout_name, rout_name=self.rout_name)
        print("Running command: {}").format(cmd)
        subprocess.check_call(cmd, shell=True)

    def test_outfile_identifiers(self):
        r1 = FastqParse(self.fout_name)
        r2 = FastqParse(self.rout_name)
        self.assertEqual(r1[self.query_fid1], True)
        self.assertEqual(r2[self.query_rid1], True)

    def test_output_num_recs(self):
        r1 = FastqParse(self.fout_name)
        r2 = FastqParse(self.rout_name)
        self.assertEqual(len(r1), 1, "Incorrect number of records in forward reads output file.")
        self.assertEqual(len(r2), 1, "Incorrect number of records in reverse reads output file.")


if __name__ == '__main__':
    unittest.main()
