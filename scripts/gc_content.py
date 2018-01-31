#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import reverse_complement, complement
from Bio.SeqUtils import GC
import numpy as np
import pyBigWig as pybw

def get_gc_coverage(fasta, binwidth, output):
    half_width = (binwidth-1)//2

    fasta_dict = {record.id:record.seq for record in SeqIO.parse(fasta, "fasta")}
    header = [(chrom, len(seq)) for chrom, seq in fasta_dict.items()]

    outbw = pybw.open(output, 'w')

    outbw.addHeader(header)
    for chrom, sequence in fasta_dict.items():
        #pad the ends of the chromosome with mirrored sequence
        padded = complement(reverse_complement(sequence[0:half_width])) + sequence + complement(reverse_complement(sequence[-half_width::]))
        #should really vectorize this
        gc_vector = np.zeros(len(sequence))
        for i in range(len(sequence)):
            gc_vector[i] = GC(padded[i:i+binwidth])
        outbw.addEntries(chrom, 0, values=gc_vector, span=1, step=1)

    outbw.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Given fasta, return bigwig of GC content using sliding window.')
    parser.add_argument('-f', dest = 'fasta', type=str, help='input fasta file')
    parser.add_argument('-w', dest = 'binwidth', type=int, help='width of sliding window for GC% calculation. MUST BE ODD INTEGER.')
    parser.add_argument('-o', dest = 'output', type=str, help='output GC% bigwig')
    args = parser.parse_args()

    get_gc_coverage(fasta=args.fasta, binwidth=args.binwidth, output=args.output)
