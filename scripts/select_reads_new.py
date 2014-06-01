#!/usr/bin/env python
import sys


def get_dictionary(in_reads):
    #erstellt dictionary
    #key = readname, value = bitflag
    d = {}

    for line in in_reads:
        read_array = str.split(line, '\t')
        if line[0] != '@':
            d[str(read_array[0])] = str(read_array[1])

    return d


def write_to_file(d, in_fq, out_reads):
    # durchlaeuft die input-fastq-datei
    # falls read in d und bitflag=4 schreiben der zum read ghoerigen zeilen in
    # die out_reads.fastq

    for line in in_fq:
        if line[0] == '@':
            id_arr = str.split(line, ' ')
            id_name_at = id_arr[0]
            id_name = id_name_at[1:]

            for paar in d.iteritems():
                if paar[0] == id_name and paar[1] == '4':
                        out_reads.write(line)
                        for i in range(3):
                            l_count = in_fq.next()
                            out_reads.write(l_count)


if __name__ == "__main__":

    infile_sam = sys.argv[1]
    infile_fq = sys.argv[2]
    outfile = sys.argv[3]

    in_reads = open(infile_sam, 'r')
    in_fq = open(infile_fq, 'r')
    out_reads = open(outfile, 'w')

    d = get_dictionary(in_reads)
    write_to_file(d, in_fq, out_reads)

    out_reads.close()
    in_reads.close()
    in_fq.close()
