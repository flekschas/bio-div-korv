# Biodiversity and Evolution - Identification of KoRV insertion sites in Koala genome (exercise 1) #

With this pipeline the KoRV specific insertion sites in the Koala genome can be identified.

### Pipeline overview ###

1. Read selection: Select reads containing the bait sequence.
2. Quality control: Trim bases with low quality and cut adapter sequences.
3. Read preperation: Discard reads that are going into the KoRV genome rather than the Koala genome and finally cut bait and remaining KoRV sequences.
4. Clustering: Cluster reads to identify different insertion sites in healthy and tumor tissue.

### Requirements ###

* Linux or Mac
* Python (>= 2.7)
* Perl
* [fqgrep](https://github.com/indraniel/fqgrep)
* [trimmomatic (>= 0.32)](http://www.usadellab.org/cms/?page=trimmomatic)
* [BWA (>= 0.7.8)](http://bio-bwa.sourceforge.net/)
* [Samtools (>= 0.1.19)](http://samtools.sourceforge.net/)
* [USearch (>= 7)](http://www.drive5.com/usearch/)

### Install and run ###

1. Clone repository `git clone https://bitbucket.org/flekschas/biodivex1 bioDivEx1`
2. Go to the project folder `cd bioDivEx1`
3. Open `scripts/run.sh` and adjust all path variables to match your system
4. Make scripts executable `chmod +x scripts/run.sh` and run `scripts/run.sh` 