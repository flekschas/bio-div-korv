#!/bin/bash

# Copyright 2014 Bressin, Kersten, Lekschas and Liedtke
# 2014-06-03
# Identify unique insertion sites for KoRV in the Koala genome.

# Requirements:
# Linux or Mac
# Python >= 2.7
# Perl
# fqgrep
# trimmomatic
# BWA
# DNAClust

# Folder structure:
#
# [BASE]
# |-data
#   |-<GROUPS>
#     |-tmp
# |-log
# |-results
#   |- 3_prime
#   |- 5_prime
# |-scripts
#
# Folder description:
#
# /data
# Contains the genome of KoRV (or any other retro virus).
# /data/<GROUPS>
# The <GROUPS> folder should contain all fastq files that belong to a group, for
# which the insertion sites are clustered. (E.g. healthy and diseased).
# /data/<GROUPS>/tmp
# The tmp folder holds all working files which can be deleted after the run.
#
# /log
# Contains log files created by trimmomatic
#
# /results
# This folder contains the cluster description given by DNAClust
# /results/3_prime
# Contains the fasta files for each cluster of the <GROUPS> found at the 3 prime
# end.
# /results/5_prime
# Contains the fasta files for each cluster of the <GROUPS> found at the 5 prime
# end.
#
# /scripts
# This folder contains all relevant scripts needed for the analysis.

################################################################################
# Config                                                                       #
################################################################################

# Please define edit all paths to match your system.
# Paths need NOT contain a trailing slash.

# Base path to your working directory
# Do NOT put a trailing slash '/' at the end of the path
BASE="/Users/Fritz/Documents/Studium/BioInformatik/Master/Module/Biodiversity and Evolution/Exercise/tool"

# Path to commands
FQGREP="fqgrep"
TRIMMOMATIC_BASE="/Applications/trimmomatic-0.32"
TRIMMOMATIC="trimmomatic-0.32.jar"
BWA="bwa"
DNACLUST="dnaclust"

# Bait sequence matching errors
BAIT_3_PRIME="CTACCCGAACATTGGGGTC"
BAIT_5_PRIME="AGGAGGCAGAAATCATGAGG"

# KoRV genome
KORV="korv"
# How many bases is the bait sequence away from the end of the genome?
KORV_3_PRIME_EXTRA=5
KORV_5_PRIME_EXTRA=5
KORV_3_PRIME_END='CTACCCGAACATTGGGGTCTTTCAT'
KORV_5_PRIME_END='AATGAAGGAGGCAGAAATCATGAGG'

# Use GNU parallel to halve running time of read selection
# At least two cores are needed to make use of parallelising
# Path to GNU parallel
PARALLEL="parallel"
# Set true to use GNU parallel
PARALLEL_USE=true



################################################################################
# Application parameters                                                       #
################################################################################

# How many mismatches should be allowed for when looking for reads that
# contain the bait sequence?
BAIT_ERR=1

# Insertion site length
INS_LEN=40

# Clustering identity threshold for DNAClust
ID=0.8



################################################################################
# Test config                                                                  #
################################################################################

# Run a few tests to avoid run time errors.

# Does base paths exists?
if [ ! -d "$BASE" ]; then
    echo "The given base path doesn't lead to a directory. Aborting!"
    echo "Path: $BASE"
    exit 1
fi

# Is the base directory writable?
if [ ! -w "$BASE" ]; then
    echo "The given base directoy exists but is not writable. Aborting!"
    echo "Path: $BASE"
    exit 1
fi

# Do applications / commands exists?
if ! type $FQGREP >/dev/null 2>&1; then
    echo "Fqgrep not found. Aborting!"
    echo "fqgrep: $FQGREP"
    exit 1
fi
# Does base paths exists?
if [ -d "$TRIMMOMATIC_BASE" ];
then
    if [ ! -f "$TRIMMOMATIC_BASE/$TRIMMOMATIC" ]; then
        echo "Trimmomatic jar file not found. Aborting!"
        echo "Trimmomatic: $TRIMMOMATIC_BASE/$TRIMMOMATIC"
        exit 1
    fi
else
    echo "The given path to trimmomatic doesn't exist. Aborting!"
    echo "Path: $TRIMMOMATIC_BASE"
    exit 1
fi
if ! type $BWA >/dev/null 2>&1; then
    echo "BWA not found. Aborting!"
    echo "BWA: $BWA"
    exit 1
fi
if ! type $DNACLUST >/dev/null 2>&1; then
    echo "DNAClust not found. Aborting!"
    echo "DNAClust: $DNACLUST"
    exit 1
fi
if $PARALLEL_USE && ! type $PARALLEL >/dev/null 2>&1; then
    echo "GNU Parallel not found. Please edit the path or install it."
    echo "Running analysis without GNU parallel."
    PARALLEL_USE=false
fi

if [ ! -f "$BASE/data/$KORV.fa" ]; then
    echo "FASTA file of the KoRV genome not found. Aborting!"
    echo "Path: $BASE/data/$KORV.fa"
    exit 1
fi

# Test if python exists
if ! type python >/dev/null 2>&1; then
    echo "Python not found. Aborting!"
    exit 1
fi

# Test if perl exists
if ! type perl >/dev/null 2>&1; then
    echo "Perl not found. Aborting!"
    exit 1
fi



################################################################################
# 1. Selection of reads containing the bait sequence                           #
################################################################################

# Start timer

ANALYSIS_START=`date +%s`

# Remove old statistics file if exists
rm -f "$BASE/results/stats.csv"

# Using fqgrep we select all reads with the defined number of errors that
# contain the two bait sequences related to the 5' and 3' LTR sequence.

echo "1. Selection of reads containing the bait sequence"

for GROUP in "$BASE/data/"*;
do
    if [ -d "${GROUP}" ]; then

        ########################################################################
        # 1. Selection of reads containing the bait sequence                   #
        ########################################################################

        # Extract the group name from the specific group path.
        GROUPNAME=$(echo ${GROUP} | perl -nle 'm/([^\/]+)$/; print $1')

        echo "Group: $GROUPNAME"

        # Create temporary folder to separate working files from original
        # source files.
        mkdir -p "${GROUP}/tmp"

        # Iterate over all fastq files
        for READS in "${GROUP}/"*".fastq"*;
        do
            # Extract file name
            # Get the string that follows the last slash ("/") until ".fastq"
            # Basically this matches [FILENAME] of a path string like so:
            # /my/path/to/my/file/[FILENAME].fastq.gz
            # or
            # /my/path/to/my/file/[FILENAME].fastq
            FILE=$(echo "${READS}" | perl -nle 'm/([^\/]+)\.fastq(\.gz)?$/; print $1')
            printf "   $FILE... "

            start=`date +%s`

            # Get reads containing the 3' and 5' LTR bait sequence
            # Skip the read selection if the file is already found!
            if $PARALLEL_USE;
            then
                # Use GNU parallel to
                [[ -f "${GROUP}/tmp/$FILE.3_prime_bait.fastq" && -f "${GROUP}/tmp/$FILE.5_prime_bait.fastq" ]] || \
                printf "$BAIT_3_PRIME\t$BAIT_ERR\t${GROUP}/tmp/$FILE.3_prime_bait.fastq\t${READS}\n$BAIT_5_PRIME\t$BAIT_ERR\t${GROUP}/tmp/$FILE.5_prime_bait.fastq\t${READS}\n" | \
                $PARALLEL --no-notice --colsep "\t" -j 2 "$FQGREP -p {1} -m {2} -o {3} {4}"
            else
                [ -f "${GROUP}/tmp/$FILE.3_prime_bait.fastq" ] || \
                $FQGREP -p $BAIT_3_PRIME -m $BAIT_ERR \
                -o "${GROUP}/tmp/$FILE.3_prime_bait.fastq" "${READS}"

                [ -f "${GROUP}/tmp/$FILE.5_prime_bait.fastq" ] || \
                $FQGREP -p $BAIT_5_PRIME -m $BAIT_ERR \
                -o "${GROUP}/tmp/$FILE.5_prime_bait.fastq" "${READS}"
            fi

            end=`date +%s`
            runtime=$((end-start))

            echo "extracted. ($runtime seconds)"
        done

        # Merge all reads for each group to one fastq file
        printf "   Merging... "

        # Remove file to avoid loop because when scripts is being run multiple times
        rm -f "${GROUP}/tmp/all.3_prime_bait.fastq"
        rm -f "${GROUP}/tmp/all.5_prime_bait.fastq"

        # Concatenate extracted reads to one fastq file
        cat "${GROUP}/tmp/"*.3_prime_bait.fastq > "${GROUP}/tmp/all.3_prime_bait.fastq"
        cat "${GROUP}/tmp/"*.5_prime_bait.fastq > "${GROUP}/tmp/all.5_prime_bait.fastq"

        echo "done!"



        ########################################################################
        # 2. Quality control                                                   #
        ########################################################################

        echo "2. Quality control"

        # Create log folder
        mkdir -p "$BASE/log"

        # Trim Illumina Adapter sequences TruSeq V3 Single Reads
        java -jar "$TRIMMOMATIC_BASE/$TRIMMOMATIC" SE -phred33 \
        -trimlog "$BASE/log/$GROUPNAME.3_prime_bait.trimming.log" \
        "${GROUP}/tmp/all.3_prime_bait.fastq" \
        "${GROUP}/tmp/all.3_prime_bait.trimmed.fastq" \
        ILLUMINACLIP:$TRIMMOMATIC_BASE/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:14 MINLEN:30 \

        java -jar "$TRIMMOMATIC_BASE/$TRIMMOMATIC" SE -phred33 \
        -trimlog "$BASE/log/$GROUPNAME.5_prime_bait.trimming.log" \
        "${GROUP}/tmp/all.5_prime_bait.fastq" \
        "${GROUP}/tmp/all.5_prime_bait.trimmed.fastq" \
        ILLUMINACLIP:$TRIMMOMATIC_BASE/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:14 MINLEN:30 \

        echo "   Done!"



        ########################################################################
        # 3. Removal of reads going into KorV                                  #
        ########################################################################

        echo "3. Removal of reads going into KorV"

        # Next we need to find all reads that align well to the KoRV genome as these
        # are the reads we are actually NOT interested in.

        printf "   Align reads... "

        # Create BWA index for the KoRV genome
        $BWA index -p "$BASE/data/$KORV" "$BASE/data/$KORV.fa"

        # Create alignment for out reads containing the bait sequences
        $BWA aln "$BASE/data/$KORV" "${GROUP}/tmp/all.3_prime_bait.trimmed.fastq" \
        > "${GROUP}/tmp/all.3_prime_bait.bwa"
        $BWA samse "$BASE/data/$KORV" "${GROUP}/tmp/all.3_prime_bait.bwa" \
        "${GROUP}/tmp/all.3_prime_bait.trimmed.fastq" \
        > "${GROUP}/tmp/all.3_prime_bait.sam"

        $BWA aln "$BASE/data/$KORV" "${GROUP}/tmp/all.5_prime_bait.trimmed.fastq" \
        > "${GROUP}/tmp/all.5_prime_bait.bwa"
        $BWA samse "$BASE/data/$KORV" "${GROUP}/tmp/all.5_prime_bait.bwa" \
        "${GROUP}/tmp/all.5_prime_bait.trimmed.fastq" \
        > "${GROUP}/tmp/all.5_prime_bait.sam"

        echo "done!"

        # Now we collect all reads that haven't matched with KoRV

        printf "   Collect unmapped reads... "

        # Create or empty a file for temporal storage of read identifiers
        >| "${GROUP}/tmp/read_ids_3_prime.exp"
        >| "${GROUP}/tmp/read_ids_5_prime.exp"

        # Get first and second column of sam file (read identifier and bitflag)
        # Search for lines not starting with an @ that did NOT matched
        # Bitflag 4 stands for unmatched reads
        # Append those reads to our temporal indentifier storage file

        cut -f 1,2 "${GROUP}/tmp/all.3_prime_bait.sam" | egrep '^[^@]' | awk '{ if($2 == 4) { print $1 } }' >> "${GROUP}/tmp/read_ids_3_prime.exp"
        cut -f 1,2 "${GROUP}/tmp/all.5_prime_bait.sam" | egrep '^[^@]' | awk '{ if($2 == 4) { print $1 } }' >> "${GROUP}/tmp/read_ids_5_prime.exp"

        echo "done!"

        # Cut all of the KoRV sequence from the reads

        printf "   Cut bait and KoRV sequence... "

        # We use fqgrep's statistics with option -r to get the starting and
        # end position of matches for.
        # For 3' we need the end match position (column 8)
        # For 5' we need the start match position (column 7)

        # $FQGREP -r -p $KORV_3_PRIME_END -m $BAIT_ERR \
        # "${GROUP}/tmp/all.3_prime_bait.trimmed.fastq" | cut -f 1,8 \
        # | "$BASE/scripts/trim_reads.py" "${GROUP}/tmp/all.3_prime_bait.trimmed.fastq" "${GROUP}/tmp/read_ids_3_prime.exp" 0 3 "${GROUP}/tmp/all.3_prime_bait_in_koala.trimmed.fasta" "fasta" $INS_LEN

        # $FQGREP -r -p $KORV_5_PRIME_END -m $BAIT_ERR \
        # "${GROUP}/tmp/all.5_prime_bait.trimmed.fastq" | cut -f 1,7 \
        # | "$BASE/scripts/trim_reads.py" "${GROUP}/tmp/all.5_prime_bait.trimmed.fastq" "${GROUP}/tmp/read_ids_5_prime.exp" 0 5 "${GROUP}/tmp/all.5_prime_bait_in_koala.trimmed.fasta" "fasta" $INS_LEN

        $FQGREP -r -p $BAIT_3_PRIME -m $BAIT_ERR \
        "${GROUP}/tmp/all.3_prime_bait.trimmed.fastq" | cut -f 1,8 \
        | "$BASE/scripts/trim_reads.py" "${GROUP}/tmp/all.3_prime_bait.trimmed.fastq" "${GROUP}/tmp/read_ids_3_prime.exp" $KORV_3_PRIME_EXTRA 3 "${GROUP}/tmp/all.3_prime_bait_in_koala.trimmed.fasta" "fasta" $INS_LEN

        $FQGREP -r -p $BAIT_5_PRIME -m $BAIT_ERR \
        "${GROUP}/tmp/all.5_prime_bait.trimmed.fastq" | cut -f 1,7 \
        | "$BASE/scripts/trim_reads.py" "${GROUP}/tmp/all.5_prime_bait.trimmed.fastq" "${GROUP}/tmp/read_ids_5_prime.exp" $KORV_5_PRIME_EXTRA 5 "${GROUP}/tmp/all.5_prime_bait_in_koala.trimmed.fasta" "fasta" $INS_LEN

        echo "done!"



        ########################################################################
        # 4. Cluster reads and create consensus sequence                       #
        ########################################################################

        # Cluster reads

        printf "4. Cluster reads..."

        # Create directory called results if it does not exist
        mkdir -p "$BASE/results"

        # Cluster reads with DNAClust using the given identity $ID
        $DNACLUST "${GROUP}/tmp/all.3_prime_bait_in_koala.trimmed.fasta" -s $ID > "$BASE/results/$GROUPNAME.3_prime.cluster.txt"
        $DNACLUST "${GROUP}/tmp/all.5_prime_bait_in_koala.trimmed.fasta" -s $ID > "$BASE/results/$GROUPNAME.5_prime.cluster.txt"

        # Create subfolders for storing clusters fasta files
        mkdir -p "$BASE/results/3_prime"
        mkdir -p "$BASE/results/5_prime"

        # Generate a fasta file for every cluster which contains all clustered
        # reads.
        "$BASE/scripts/cluster_to_fasta.py" "$BASE/results/$GROUPNAME.3_prime.cluster.txt" "${GROUP}/tmp/all.3_prime_bait_in_koala.trimmed.fasta" "$BASE/results" "3_prime" $GROUPNAME
        "$BASE/scripts/cluster_to_fasta.py" "$BASE/results/$GROUPNAME.5_prime.cluster.txt" "${GROUP}/tmp/all.5_prime_bait_in_koala.trimmed.fasta" "$BASE/results" "5_prime" $GROUPNAME

        echo "done!"
    fi
done

ANALYSIS_END=`date +%s`
ANALYSIS_RUNTIME=$((ANALYSIS_END-ANALYSIS_START))

echo "Finished! ($ANALYSIS_RUNTIME seconds)"
