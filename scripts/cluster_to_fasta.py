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
    faDir = sys.argv[3]  # output folder
    prefix = sys.argv[4]  # prefix

    cluster = open(clusterFile, "r")
    origin_fasta = open(originFastaFile, "r")
    stat_file = open(faDir + "stat.txt", 'w')

    stat_file.write("@ ClusterName \t SequencesPerCluster \n")
    dictionary = get_dictionary(origin_fasta)

    num = 0
    sequence = ""
    seq_count = 0
    for line in cluster:
        clusterA = line.strip().split('\t')
        stat_file.write(prefix + str(num) + "\t" +
                        str(len(clusterA)) + "\n")

        with open(faDir + prefix + "." + str(num) + ".fasta", "w") as fq:
            seq_count += len(clusterA)
            for i in clusterA:
                if i in dictionary:
                    header = ">" + i + "\n"
                    sequence = dictionary[i] + "\n"
                    fq.write(header)
                    fq.write(sequence)
        num += 1

    stat_file.write("found " + str(num) + " cluster with " +
                    str(seq_count) + " sequences")

    cluster.close()
    origin_fasta.close()
    stat_file.close()
