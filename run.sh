#!/bin/bash

# Fritz Lekschas,
# 2014-05-xx
# Identify unique insertion sites for KoRV in the Koala genome.

# Requirements:
# BWA
# fqgrep
# trimmomatic

# Folder structure:
# [BASE]
# |-data
#   |-<GROUPS>
#     |-tmp
# |-results
# |-scripts
#
# The <GROUPS> folder should contain all fastq files that belong to a group for
# which the insertion sites should be clustered. (E.g. healthy and diseased)

################################################################################
# Config                                                                       #
################################################################################

# Please define edit all paths to match your system.
# Paths need NOT contain a trailing slash.

# Base path to your working directory
BASE="/Users/Fritz/Documents/Studium/BioInformatik/Master/Module/Biodiversity and Evolution/Exercise"

# Path to commands
TRIMMOMATIC="/Applications/trimmomatic-0.32"
VELVET="/Applications/velvet"
USEARCH="/Applications/usearch"

# Bait sequence matching errors
BAIT_3_PRIME="CTACCCGAACATTGGGGTC"
BAIT_5_PRIME="AGGAGGCAGAAATCATGAGG"

# KoRV genome
KORV="korv"
KORV_3_PRIME_END='CTACCCGAACATTGGGGTCTTTCAT'
KORV_5_PRIME_END='AATGAAGGAGGCAGAAATCATGAGG'

# Use GNU parallel to massively speed up read selection
# At least two cores are needed to make use of parallelising the selection
PARALLEL=true



################################################################################
# Application parameters                                                       #
################################################################################

# How many mismatches should be allowed for when looking for reads that
# contain the bait sequence?
BAIT_ERR=1

# Insertion site length
INS_LEN=40

# Clustering identity threshold
ID=0.8


################################################################################
# 0 . Test config                                                              #
################################################################################

# Run a few tests:
# paths exists?
# application exists?
# python?
# perl?



################################################################################
# 1. Selection of reads containing the bait sequence                           #
################################################################################

# Using fqgrep we select all reads with the defined number of errors that
# contain the two bait sequences related to the 5' and 3' LTR sequence.

echo "1. Selection of reads containing the bait sequence"

for GROUP in "$BASE/data/"*;
do
    if [ -d "${GROUP}" ]; then
        echo "${GROUP}" | perl -nle 'm/([^\/]+)$/; print "   Group: $1"'
        # Create temporary folder to separate working files from original
        # source files.
        mkdir -p "${GROUP}/tmp"


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
            if $PARALLEL;
            then
                # Use GNU parallel to
                [[ -f "${GROUP}/tmp/$FILE.3_prime_bait.fastq" && -f "${GROUP}/tmp/$FILE.5_prime_bait.fastq" ]] || \
                printf "$BAIT_3_PRIME\t$BAIT_ERR\t${GROUP}/tmp/$FILE.3_prime_bait.fastq\t${READS}\n$BAIT_5_PRIME\t$BAIT_ERR\t${GROUP}/tmp/$FILE.5_prime_bait.fastq\t${READS}\n" | \
                parallel --no-notice --colsep "\t" -j 2 "fqgrep -p {1} -m {2} -o {3} {4}"
            else
                [ -f "${GROUP}/tmp/$FILE.3_prime_bait.fastq" ] || \
                fqgrep -p $BAIT_3_PRIME -m $BAIT_ERR \
                -o "${GROUP}/tmp/$FILE.3_prime_bait.fastq" "${READS}"

                [ -f "${GROUP}/tmp/$FILE.5_prime_bait.fastq" ] || \
                fqgrep -p $BAIT_5_PRIME -m $BAIT_ERR \
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
        cat "${GROUP}/tmp/"*.3_prime_bait.fastq > "${GROUP}/tmp/all.3_prime_bait.fastq"
        rm -f "${GROUP}/tmp/all.5_prime_bait.fastq"
        cat "${GROUP}/tmp/"*.5_prime_bait.fastq > "${GROUP}/tmp/all.5_prime_bait.fastq"

        echo "done!"

    fi
done

exit

# Get reads containing the 3' and 5' LTR bait sequence
fqgrep -p $BAIT_3_PRIME -m $BAIT_ERR \
-o "$BASE/data/$TEST_FILE.3_prime_bait.fastq" "$BASE/data/$TEST_FILE.fastq.gz"

fqgrep -p $BAIT_5_PRIME -m $BAIT_ERR \
-o "$BASE/data/$TEST_FILE.5_prime_bait.fastq" "$BASE/data/$TEST_FILE.fastq.gz"

################################################################################
# 2. Quality control                                                           #
################################################################################

echo "2. Quality control"

# Trim Illumina Adapter sequences TruSeq V3 Single Reads
java -jar $TRIMMOMATIC/trimmomatic-0.32.jar SE -phred33 \
-trimlog "$BASE/log/$TEST_FILE.3_prime_bait.trimming.log" \
"$BASE/data/$TEST_FILE.3_prime_bait.fastq" \
"$BASE/data/$TEST_FILE.3_prime_bait.trimmed.fastq" \
ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:14 MINLEN:30

java -jar $TRIMMOMATIC/trimmomatic-0.32.jar SE -phred33 \
-trimlog "$BASE/log/$TEST_FILE.5_prime_bait.trimming.log" \
"$BASE/data/$TEST_FILE.5_prime_bait.fastq" \
"$BASE/data/$TEST_FILE.5_prime_bait.trimmed.fastq" \
ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:14 MINLEN:30

# Filter reads again in case the bait sequence was found by accident but the bait
# Sequence is represented by a too low score [optional]



################################################################################
# 3. Removal of reads going into KorV                                          #
################################################################################

echo "3. Removal of reads going into KorV"

# Next we need to find all reads that align well to the KoRV genome as these
# are the reads we are actually NOT interested in.

# Create BWA index for the KoRV genome
bwa index -p "$BASE/data/$KORV" "$BASE/data/$KORV.fa"

# Create alignment for out reads containing the bait sequences
bwa aln "$BASE/data/$KORV" "$BASE/data/$TEST_FILE.3_prime_bait.fastq" \
> "$BASE/data/$TEST_FILE.3_prime_bait.bwa"
bwa samse "$BASE/data/$KORV" "$BASE/data/$TEST_FILE.3_prime_bait.bwa" \
"$BASE/data/$TEST_FILE.3_prime_bait.fastq" \
> "$BASE/data/$TEST_FILE.3_prime_bait.sam"

bwa aln "$BASE/data/$KORV" "$BASE/data/$TEST_FILE.5_prime_bait.fastq" \
> "$BASE/data/$TEST_FILE.5_prime_bait.bwa"
bwa samse "$BASE/data/$KORV" "$BASE/data/$TEST_FILE.5_prime_bait.bwa" \
"$BASE/data/$TEST_FILE.5_prime_bait.fastq" \
> "$BASE/data/$TEST_FILE.5_prime_bait.sam"

# Now we need to remove the reads contained in the sam file from the
# originally extracted reads that contained the bait sequence.

# Create or empty a file for temporal storage of reads identifier
>| "$BASE/data/tmp_read_ids_3_prime.exp"
>| "$BASE/data/tmp_read_ids_5_prime.exp"

# Get first and second column of sam file (read identifier and bitflag)
# Search for lines not starting with an @ that are NOT unmatched
# Bitflag 4 stands for unmatched reads
# Append those reads to our temporal indentifier storage file
cut -f 1,2 "$BASE/data/$TEST_FILE.3_prime_bait.sam" | egrep '^[^@]' | awk '{ if($2 != 4) { print "/" $1 "/{N;N;N;d;}" } }' >> "$BASE/data/tmp_read_ids_3_prime.exp"
cut -f 1,2 "$BASE/data/$TEST_FILE.5_prime_bait.sam" | egrep '^[^@]' | awk '{ if($2 != 4) { print "/" $1 "/{N;N;N;d;}" } }' >> "$BASE/data/tmp_read_ids_5_prime.exp"

# Delete all reads that are not unmatched using sed
# sed loads the temporal read storage file. Each line contains a string as follows
# Syntax: /<READ-ID>/{N;N;d;} matches lines starting with <READ-ID> plus another
# 3 lines after a match. The d flag says sed should delete those lines.
sed -f "$BASE/data/tmp_read_ids_3_prime.exp" "$BASE/data/$TEST_FILE.3_prime_bait.fastq" > "$BASE/data/$TEST_FILE.3_prime_bait_no_korv.fastq"
sed -f "$BASE/data/tmp_read_ids_5_prime.exp" "$BASE/data/$TEST_FILE.5_prime_bait.fastq" > "$BASE/data/$TEST_FILE.5_prime_bait_no_korv.fastq"

# # Prepare FASTQ files for contig assembly
# $VELVET/velveth "$BASE/data/velvet_3_prime/" 19 -short -fastq "$BASE/data/$TEST_FILE.3_prime_bait.fastq"
# $VELVET/velveth "$BASE/data/velvet_5_prime/" 19 -short -fastq "$BASE/data/$TEST_FILE.5_prime_bait.fastq"

# # Assemble reads
# $VELVET/velvetg "$BASE/data/velvet_3_prime" -exp_cov auto
# $VELVET/velvetg "$BASE/data/velvet_5_prime" -exp_cov auto

# Cut all of the KoRV sequence from the reads
# -r prints match statistics
# For 3' we need the end match position (column 8)
# For 5' we need the start match position (column 7)
fqgrep -r -p $KORV_3_PRIME_END -m $BAIT_ERR \
"$BASE/data/$TEST_FILE.3_prime_bait_no_korv.fastq" | cut -f 1,8 \
| "$BASE/scripts/trim_reads.py" "$BASE/data/$TEST_FILE.3_prime_bait_no_korv.fastq" 3 "$BASE/data/$TEST_FILE.3_prime_bait_no_korv_trimmed.fasta" "fasta" $INS_LEN

fqgrep -r -p $KORV_5_PRIME_END -m $BAIT_ERR \
"$BASE/data/$TEST_FILE.5_prime_bait_no_korv.fastq" | cut -f 1,7 \
| "$BASE/scripts/trim_reads.py" "$BASE/data/$TEST_FILE.5_prime_bait_no_korv.fastq" 5 "$BASE/data/$TEST_FILE.5_prime_bait_no_korv_trimmed.fasta" "fasta" $INS_LEN

# Cluster reads
# Create directory called results if it does not exist
mkdir -p "$BASE/results"

# Clustering
$USEARCH -cluster_fast "$BASE/data/$TEST_FILE.3_prime_bait_no_korv_trimmed.fasta" -id $ID \
-centroids "$BASE/results/$TEST_FILE.3_prime_"$INS_LEN"BP_cluster.fasta"
$USEARCH -cluster_fast "$BASE/data/$TEST_FILE.5_prime_bait_no_korv_trimmed.fasta" -id $ID \
-centroids "$BASE/results/$TEST_FILE.5_prime_"$INS_LEN"BPcluster.fasta"
