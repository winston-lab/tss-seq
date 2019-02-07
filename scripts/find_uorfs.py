#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import re

def main(fasta_path: str,
         input_bed_path: str,
         output_path: str):
    stop_codon_regex = "TAA|TAG|TGA"
    bed = pd.read_csv(input_bed_path,
                      sep="\t",
                      names=["chrom", "start", "end", "name", "score", "strand"])
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    with open(output_path, "w") as output_file:
        for utr_index, utr in bed.iterrows():
            utr_sequence = fasta[utr.chrom][utr.start : utr.end]
            if utr.strand == "-":
                utr_sequence = utr_sequence.reverse_complement()

            atg_position = utr_sequence.seq.find("ATG")
            uorf_counter = 0
            while atg_position != -1:
                stop_position = [stop_codon.end() for stop_codon in \
                        re.finditer(stop_codon_regex, str(utr_sequence.seq)) if \
                        stop_codon.start() > atg_position and \
                        (stop_codon.start() - atg_position) % 3 == 0]
                if stop_position:
                    stop_position = stop_position[0]
                    uorf_counter += 1
                    if utr.strand == "+":
                        uorf_start = utr.start + atg_position
                        uorf_end = utr.start + stop_position
                    elif utr.strand == "-":
                        uorf_start = utr.end - stop_position
                        uorf_end = utr.end - atg_position
                    output_file.write("\t".join([utr.chrom,
                                                 str(uorf_start),
                                                 str(uorf_end),
                                                 utr["name"] + "_uORF_" + str(uorf_counter),
                                                 str(0),
                                                 utr.strand]) + "\n")
                atg_position = utr_sequence.seq.find("ATG", atg_position+1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given a BED6 file of 5\'-UTRs and a FASTA of the genome sequence, find uORFs.')
    parser.add_argument('-f', dest = 'fasta_path', type=str, help='FASTA file of genome.')
    parser.add_argument('-i', dest = 'input_bed_path', type=str, help='BED6 file of 5\'UTR coordinates.')
    parser.add_argument('-o', dest = 'output_path', type=str, help='Output BED6 file uORF coordinates.')
    args = parser.parse_args()

    main(fasta_path=args.fasta_path,
         input_bed_path=args.input_bed_path,
         output_path=args.output_path)

