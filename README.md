# crisprcas9_gRNA_design
My solution to the PhD position assignment in the Gorodkin lab.

This script can look for 5'UTR in a given FASTA file using a GTF annotation file, and design guide RNAs to target these regions.

This script can be run from the terminal. There are 3 required arguments:
 1. Fasta file containing the reference genome
 2. GTF file containing genome annotation
 3. Number of the chromosome to work with

## Example command for running the script:
> python design_5utr_guideRNA.py GRCh38.primary_assembly.genome.fa gencode.v35.annotation.gtf 22
