# illumina_fastq

Contains a module for parsing a FASTQ file and several scripts based upon it. The present scripts are:

#### illuminaBarcodeDist.py
Tabulates the relative frequences of barcodes seen in the unmatched reads. This is useful when you have a lot of unmatched reads and need to see what is present. 

There is a Docker image for this at https://hub.docker.com/r/nathankw/illumina_barcode_dist. 

#### extract_barcodes.py
Extracts FASTQ records matching the specified barcodes from the unmatched read FASTQ file or pair of FASTQ files if paired-end sequencing.

There is a Docker image for this at https://hub.docker.com/r/nathankw/illumina_extract_barcodes. 
#### add_mate_reads.py
Fetches the other read in a pair for each record in the query FASTQ file. You also supply the parent pair of forward read and reverse read FASTQ files as the reference set of FASTQ records.

#### get_ids_in_fastq_files.py
Generates a file containing the FASTQ record identifiers (one per line) from the passed in FASTQ files.

### Unit Tests
There are tests modeled after the scripts within the folder scripts/Test within a subfolder named after each script. Each test script can be run without arguments. 
