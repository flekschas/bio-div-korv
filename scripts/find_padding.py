#!/usr/bin/env python
import sys


def remove_lines(filename, line_numbers, num_of_lines):
    lines = []
    skip_lines = 0
    with open(line_numbers, 'r') as f:
        for line in f:
            lines.append(int(line.strip()))
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i in lines:
                skip_lines = i + num_of_lines
                print(i, skip_lines)
            elif i > skip_lines:
                print(line.strip())


#############
# MAIN      #
#############

def main():
    if len(sys.argv) > 3:
        remove_lines(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        print("rm_lines.py <fastq> <line number>")
    sys.exit()

if __name__ == "__main__":
    main()
