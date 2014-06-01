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
    #fq_file = open(fastaFile, "w")

    dictionary = get_dictionary(origin_fasta)

    number = 0
    sequence = ""
    for line in cluster:
        # line = line.strip()
        clusterA = line.strip().split('\t')
        StringName = str(fastaFile_folder) + str(outputName) + str(number) + ".fasta"
        fq_file = open(StringName, "w")
        for i in clusterA:
            if i in dictionary:
                header = ">" + i + "\n"
                sequence = dictionary[i] + "\n"
                fq_file.write(header)
                fq_file.write(sequence)
            else:
                print('hallo')
        number += 1

    cluster.close()
    origin_fasta.close()
    fq_file.close()
