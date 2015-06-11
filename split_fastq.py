#!/usr/bin/env python
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Split interleaved FastQ')
parser.add_argument('--trim', nargs=1, type=int, default=0,
                    help='Trim N chars from trailing end of sequence names')
parser.add_argument('input', metavar='FASTQ', help='Input FASTQ format sequence to split')
parser.add_argument('output', metavar='BASE_OUT', help='Basename for output files')
args = parser.parse_args()

r1_out = '{0}1.fq'.format(args.output[0])
r2_out = '{0}2.fq'.format(args.output[0])

with open(r1_out, 'w') as r1_h, open(r2_out, 'w') as r2_h:

    seq_it = SeqIO.parse(args.input[0], 'fastq')
    c1 = c2 = 0

    try:
        if args.trim > 0:
            nch = args.trim
            while True:
                seq = seq_it.next()
                seq.id = seq.id[:-nch]
                SeqIO.write([seq], r1_h, 'fastq')
                c1 += 1
                seq = seq_it.next()
                seq.id = seq.id[:-nch]
                SeqIO.write([seq], r2_h, 'fastq')
                c2 += 1

        else:
            while True:
                SeqIO.write([seq_it.next()], r1_h, 'fastq')
                c1 += 1
                SeqIO.write([seq_it.next()], r2_h, 'fastq')
                c2 += 1

    except StopIteration:
        print 'Wrote {0} forward and {1} reverse reads'.format(c1, c2)
        if c1 != c2:
            print 'Warning: forward and reverse sets were not of equal length'
