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
* Perl (>= 5.12)
* [fqgrep](https://github.com/indraniel/fqgrep)
* [trimmomatic (>= 0.32)](http://www.usadellab.org/cms/?page=trimmomatic)
* [BWA (>= 0.7.8)](http://bio-bwa.sourceforge.net/)
* [DNAClust (>= 3)](http://dnaclust.sourceforge.net/)
* [GNU Parallel](https://www.gnu.org/software/parallel/) (optional)

### Install and run ###

1. Clone repository `git clone https://bitbucket.org/flekschas/biodivex1 bioDivEx1`
2. Go to the project folder `cd bioDivEx1`
3. Open `nano scripts/run.bash` and edit all path variables to match your system which can be found in section *config*.
4. Adjust the application parameters to suite your problem in section *application parameters*.
5. Make scripts executable `chmod +x scripts/cluster_to_fasta.py`, `chmod +x scripts/run.bash`, `chmod +x scripts/trim_reads.py`.
6. Copy your `FASTQ` files into the data folder and group them into subfolders for later comparison. E.g. all `FASTQ` files from healthy tissues could be placed into `data/healthy` and all FASTQ files related to diseased tissues respecively into `data/diseased`. For detailed
explanations please have a look at the `README.md` files inside the data folder.
7. Copy the `FASTA` file with your retro virus genome into `data`. For the KoRV analysis we already included `korv.fa`.
8. Run `scripts/runs.bash`

### Results ###

When the execution finished your should see a new folder called `results`. Within this folder you can find a clustering file of your possible insertion sites for 3 prime and 5 prime reads of all of your groups. E.g. [GROUPNAME].3_prime.cluster.txt and [GROUPNAME].5_prime.cluster.txt.

Two subfolders called `results/3_prime` and `results/5_prime` hold `FASTA` files for each of the found instertion site clusters.

Finally you should find `results/stats.csv` a small file holding a few statistics about the clusters.

### Run time ###

The first run takes about 15 minutes on a MacBook Air 2012 with 2 cores and GNU parallel enabled. The bait selection
takes most of the time.
A second run on the same data using the same groups then only takes about 6 seconds as the pipeline looks for already
selected reads. So parameter adjustments is confenient.

**Note**: If you change the bait sequence and still want to use the old sequence files and groups please delete the `data/[GROUP]/tmp` folders first. Otherwise the script won't use the new bait sequence.
