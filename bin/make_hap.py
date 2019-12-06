#!/usr/bin/env python2

from random import random, choice
import sys
from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def main(fasta1, fasta2):
    iter1 = fasta_iter(fasta1)
    iter2 = fasta_iter(fasta2)
    iter1_empty = False
    iter2_empty = False
    while True:
        try:
            header1, seq1 = iter1.next()
        except StopIteration:
            iter1_empty = True
        try:
            header2, seq2 = iter2.next()
        except StopIteration:
            iter2_empty = True

        assert iter1_empty == iter2_empty, "Input FASTA files have different number of sequences."

        if iter1_empty and iter2_empty:
            break

        assert header1 == header2, "Headers don't match; both input FASTAs should be sorted and have the same headers."

        header = header1

        assert len(seq1) == len(seq2), "Length of sequences do not match."

        hap = ['N'] * len(seq1)

        for i in range(len(seq1)):
            s1 = seq1[i]
            s2 = seq2[i]
            hap[i] = choice([s1, s2])
        print(">%s\n%s" % (header, "".join(hap)))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
