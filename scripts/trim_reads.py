#!/usr/bin/env python
import sys


def trim_reads(fastq, selection, orientation, output, outputType, seqLen,
               trim):

    # Store all read/sequence ids that did not match with KoRV
    ids = []
    with open(selection, 'r') as f:
        for line in f:
            ids.append(line.strip())

    # Store trimming position for each read/sequence id
    trim_pros = {}
    for line in trim.split('\n'):
        if len(line):
            line = line.split('\t')
            if (line[0] == 'read name'):
                if (line[1] == 'end position' and orientation != 3) or \
                   (line[1] == 'start position' and orientation != 5):
                    print('Wrong setting! 3\' trimming needs the end position'
                          'and 3\' trimming needs the start position.')
                    sys.exit()
            else:
                trim_pros[line[0]] = int(line[1])

    # Read fastq file line by line and copy a sequence to a new fastq file if:
    # 1. Read did not align against KoRV (id is in selection)
    # 2. Line is not blank
    # 3. Sequence length is greater than the given seqLen
    with open(output, 'w') as o:
        with open(fastq, 'r') as f:
            while True:
                identifier = f.readline()
                sequence = f.readline()
                plus = f.readline()
                quality = f.readline()

                if not identifier or not sequence or \
                   not plus or not quality:
                    break

                read_id = identifier.strip()[1:].split(' ')[0]

                if (read_id in ids and read_id in trim_pros):
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
    trim = sys.stdin.read()
    if len(sys.argv) > 6:
        trim_reads(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4],
                   sys.argv[5], int(sys.argv[6]), trim)
    else:
        print("trim_reads.py [fastq] [selection] [orientation] [output] "
              "[format] [maxlen] < [trimming-info]")
    sys.exit()

if __name__ == "__main__":
    main()
