# -*- coding: utf-8 -*-                                                                                
                                                                                                       
###                                                                                                    
# © 2018 The Board of Trustees of the Leland Stanford Junior University                                
# Nathaniel Watson                                                                                     
# nathankw@stanford.edu                                                                                
###

# For some useful documentation, see
# https://docs.python.org/2/distutils/setupscript.html.
# This page is useful for dependencies: 
# http://python-packaging.readthedocs.io/en/latest/dependencies.html.

import glob
import os
from setuptools import setup, find_packages

SCRIPTS_DIR = "illumina_fastq/scripts/"
scripts = glob.glob(os.path.join(SCRIPTS_DIR,"*.py"))
scripts.remove(os.path.join(SCRIPTS_DIR,"__init__.py"))

setup(
  name = "illumina_fastq",
  version = "1.0.0",
  description = "Tools for parsing Illumina FASTQ files. You can create a barcode distribution of unmatched reads, and extract barcodes of interest.",
  author = "Nathaniel Watson",
  author_email = "nathankw@stanford.edu",
  url = "https://github.com/StanfordBioinformatics/illumina_fastq",
  packages = find_packages(),
  install_requires = [],
  scripts = scripts,
  package_data = {"encode_utils": ["tests/data/*"]}
)
