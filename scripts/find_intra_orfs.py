#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import molecular_weight

parser = argparse.ArgumentParser(description='Detect potential ORFs downstream of intragenic TSS-seq peaks.')
parser.add_argument('-p', dest = 'peaks', type=str, default = 'diamide-v-YPD-de-clusters-libsizenorm-up-intragenic.tsv')
parser.add_argument('-f', dest = 'fasta', type=str, default = 'S_cerevisiae.R64-2-1.fa')
parser.add_argument('-m', dest = 'max_orf_len', type=int, default = 10000)
parser.add_argument('-o', dest = 'outpath', type=str, default = 'output.tsv')
args = parser.parse_args()

peaks = pd.read_table(args.peaks, names = ['chrom','strand','peakstart','peakend', 'peakname', 'orfstart', 'orfend', 'orfname', 'peaksig', 'ATGtoPeakDist'])
fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

peaks['frame'] = np.where(peaks['strand']=="+", (peaks['peakstart']-peaks['orfstart']) % 3, (peaks['orfend']-peaks['peakend']) % 3)

#initialize output dataframe
df = pd.DataFrame(index = np.arange(len(peaks)*1000), 
                  columns = ['chrom', 'strand', 'orf_name', 'orf_start', 'orf_end', 'dist_atg_to_peak', 'peak_name', 'peak_start', 'peak_end', 'peak_significance', 'intra_orf_5utr_length',  'intra_orf_frame', 'intra_orf_upstr_atg', 'intra_orf_start', 'intra_orf_end', 'intra_orf_length', 'intra_prot_molweight', 'intra_prot_tap_molweight'])

tt_out = 0 #counter for output df row

#make TAP tag sequence object
tap_tag = Seq('SMEKRRWKKNFIAVSAANRFKKISSSGALDYDIPTTASENLYFQGELKTAALAQHDEAVDNFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPKVDNFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPK',                IUPAC.extended_protein)
tap_molweight = molecular_weight(tap_tag)

for index, row in peaks.iterrows():
    atg_counter = 0

    if row.strand == "+":
        seq = fasta[str(row.chrom)].seq[row.peakstart : row.orfend + args.max_orf_len]
        
        startpos = seq.find("ATG")
        intra_orf_start = row.peakstart
        frame = row.frame
        
        while (intra_orf_start < row.orfend and startpos > -1):
            if startpos > -1:
                intra_orf_start = row.peakstart + startpos
                frame = (frame + startpos) % 3 #get reading frame of found start codon

                intraprotein = seq[startpos:].translate(to_stop=True)
                ntlen = (len(intraprotein) + 1) * 3

                if (intraprotein.find("X") == -1):
                    molecweight = molecular_weight(intraprotein)
                else:
                    molecweight = -1

                df.iloc[tt_out,] = [row.chrom, row.strand, row.orfname, row.orfstart, row.orfend, row.ATGtoPeakDist, row.peakname, row.peakstart, row.peakend, row.peaksig, startpos, frame, atg_counter, intra_orf_start, intra_orf_start+ntlen, ntlen,  molecweight/1000, (molecweight+tap_molweight)/1000]
                tt_out += 1
                atg_counter += 1
                startpos = seq.find("ATG", start = startpos + 1)
            else:
                intra_orf_start = row.orfend


    elif row.strand == "-":
        seq = fasta[str(row.chrom)].seq[row.orfstart - args.max_orf_len : row.peakend].reverse_complement()

        startpos = seq.find("ATG")
        intra_orf_start = row.peakend 
        frame = row.frame

        while (intra_orf_start > row.orfstart):
            if startpos > -1:
                intra_orf_start = row.peakend - startpos
                frame = (frame + startpos) % 3 #get reading frame of found start codon

                intraprotein = seq[startpos:].translate(to_stop=True)
                ntlen = (len(intraprotein) + 1) * 3

                if (intraprotein.find("X") == -1):
                    molecweight = molecular_weight(intraprotein)
                else:
                    molecweight = -1

                df.iloc[tt_out,] = [row.chrom, row.strand, row.orfname, row.orfstart, row.orfend, row.ATGtoPeakDist, row.peakname, row.peakstart, row.peakend, row.peaksig, startpos, frame, atg_counter, intra_orf_start-ntlen, intra_orf_start, ntlen, molecweight/1000, (molecweight+tap_molweight)/1000]
                tt_out += 1
                atg_counter += 1
                startpos = seq.find("ATG", start = startpos + 1)
            else:
                intra_orf_start = row.orfstart   

df = df.iloc[:tt_out]

df.to_csv(args.outpath, sep = "\t", header=True, index=False)
