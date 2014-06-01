#!/usr/bin/env python
import sys


def select_reads(fastq, output, selection):
    selection = selection.split('\n')

    with open(output, 'w') as o:
        with open(fastq, 'r') as f:
            while True:
                identifier = f.readline()
                sequence = f.readline()
                plus = f.readline()
                quality = f.readline()

                if (not identifier or
                    not sequence or
                    not plus or
                    not quality):
                    break

                read_id = identifier.strip()[1:].split(' ')[0]

                if (read_id in trim_pros):
                    if (orientation == 3):
                        sequence = sequence[trim_pros[read_id]:(trim_pros[read_id] + seqLen)].strip()
                        quality = quality[trim_pros[read_id]:(trim_pros[read_id] + seqLen)].strip()
                    if (orientation == 5):
                        sequence = sequence[max(trim_pros[read_id] - seqLen, 0):trim_pros[read_id]]
                        quality = quality[max(trim_pros[read_id] - seqLen, 0):trim_pros[read_id]]

                    if (len(sequence) >= seqLen):
                        if (outputType == 'fasta'):
                            o.write('>' + identifier[1:])
                            o.write(sequence[:seqLen] + '\n')
                        else:
                            o.write(identifier)
                            o.write(sequence[:seqLen] + '\n')
                            o.write(plus)
                            o.write(quality[:seqLen] + '\n')


#############
# MAIN      #
#############

def main():
    selection = sys.stdin.read()
    if len(sys.argv) > 2:
        select_reads(sys.argv[1], sys.argv[2], selection)
    else:
        print("select_reads.py [fastq] [output] < [trimming-info]")
    sys.exit()

if __name__ == "__main__":
    main()
