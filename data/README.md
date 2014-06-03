# Data #

This folder contains the `FASTA` files of our retro virus KoRV and subfolders
for each of the groups we want to compare. In our case `healthy` and `diseased`
tissue.

## KoRV: Folder structure ##

Our folder structure for analysing KoRV insertion sites looks like follows.

### Healthy vs tumour tissue ###

In this case we pool all healthy and tumour samples together and compare those
two groups only.

* data/
** diseased/
*** TOG-ECAK-47628-T_S4_L001_R1_001.fastq.gz
*** TOG-ECAK-47813-T_S3_L001_R1_001.fastq.gz
*** TOG-ECAK-A30038-T_S6_L001_R1_001.fastq.gz
*** TOG-ECAK-ELSIE-T_S5_L001_R1_001.fastq.gz
*** TOG-ECAK-MG-T_S1_L001_R1_001.fastq.gz
** healthy/
*** TOG-ECAK-47628-H_S4_L001_R1_001.fastq.gz
*** TOG-ECAK-47813-H_S3_L001_R1_001.fastq.gz
*** TOG-ECAK-A30038-H_S6_L001_R1_001.fastq.gz
*** TOG-ECAK-ELSIE-H_S5_L001_R1_001.fastq.gz
*** TOG-ECAK-MG-H_S1_L001_R1_001.fastq.gz
** korv.fa

### Sample vs sample ###

We can rearrange the folders like follows to compare all of the samples against
each other.

* data/
** TOG-ECAK-47628-T_S4_L001_R1_001/
*** TOG-ECAK-47628-T_S4_L001_R1_001.fastq.gz
** TOG-ECAK-47813-T_S3_L001_R1_001/
*** TOG-ECAK-47813-T_S3_L001_R1_001.fastq.gz
** TOG-ECAK-A30038-T_S6_L001_R1_001/
*** TOG-ECAK-A30038-T_S6_L001_R1_001.fastq.gz
** TOG-ECAK-ELSIE-T_S5_L001_R1_001/
*** TOG-ECAK-ELSIE-T_S5_L001_R1_001.fastq.gz
** TOG-ECAK-MG-T_S1_L001_R1_001/
*** TOG-ECAK-MG-T_S1_L001_R1_001.fastq.gz
** TOG-ECAK-47628-H_S4_L001_R1_001/
*** TOG-ECAK-47628-H_S4_L001_R1_001.fastq.gz
** TOG-ECAK-47813-H_S3_L001_R1_001/
*** TOG-ECAK-47813-H_S3_L001_R1_001.fastq.gz
** TOG-ECAK-A30038-H_S6_L001_R1_001/
*** TOG-ECAK-A30038-H_S6_L001_R1_001.fastq.gz
** TOG-ECAK-ELSIE-H_S5_L001_R1_001/
*** TOG-ECAK-ELSIE-H_S5_L001_R1_001.fastq.gz
** TOG-ECAK-MG-H_S1_L001_R1_001/
*** TOG-ECAK-MG-H_S1_L001_R1_001.fastq.gz
** korv.fa
