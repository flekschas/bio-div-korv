#!/usr/bin/env python
import sys


def get_dictionary(in_reads):
    #erstellt dictionary
    #key = readname, value = bitflag
    d = {}
    identify = ""
    for line in in_reads:
        line = line.strip()
        if line[0] == '>':
            identify = line[1:].split(' ')[0]
        else:
            d[identify] = line
    return d


if __name__ == "__main__":
    clusterFile = sys.argv[1]
    originFastaFile = sys.argv[2]
    fastaFile_folder = sys.argv[3]  # output folder
    outputName = sys.argv[4]  # prefix

    cluster = open(clusterFile, "r")
    origin_fasta = open(originFastaFile, "r")
    stat_file = open(fastaFile_folder + "stat.txt", 'w')
    #fq_file = open(fastaFile, "w")

    stat_file.write("@ ClusterName \t SequencesPerCluster \n" )
    dictionary = get_dictionary(origin_fasta)

    number = 0
    sequence = ""
    seq_count = 0
    for line in cluster:
        # line = line.strip()
        clusterA = line.strip().split('\t')
        StringName = str(fastaFile_folder) + str(outputName) + str(number) + ".fasta"
        fq_file = open(StringName, "w")
        stat_file.write(outputName + str(number) + "\t" + str(len(clusterA)) + "\n")
        seq_count += len(clusterA)
        for i in clusterA:
            if i in dictionary:
                header = ">" + i + "\n"
                sequence = dictionary[i] + "\n"
                fq_file.write(header)
                fq_file.write(sequence)
            else:
                print('hallo')
        number += 1
    stat_file.write("found " + str(number) +" cluster with "+ str(seq_count) + " sequences")

    cluster.close()
    origin_fasta.close()
    fq_file.close()
