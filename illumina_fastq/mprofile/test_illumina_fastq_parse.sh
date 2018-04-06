#!/bin/bash -eu 

#$ -M nathankw@stanford.edu
#$ -m ae
#$ -cwd
#$ -R y
#$ -l h_vmem=35G

module load python/2.7.9
mprof run illumina_fastq_parse.py /srv/gsfs0/projects/gbsc/workspace/nathankw/CIRM/CORN/TEST/25mil_recs.fastq
