#!/usr/bin/env python

import argparse
import pybedtools as pybt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Make BED of genic regions given transcript BED and ORF bed.')
parser.add_argument('-t', dest = 'transcripts', type=str, default = 'Scer_polIItranscripts-adjustedTSS.bed')
parser.add_argument('-o', dest = 'orfs', type=str, default = 'Scer_nondubious_ORFs.bed')
parser.add_argument('-d', dest = 'dist', type=int, default = 30)
parser.add_argument('-g', dest = 'chrsizes', type=str, default = 'S_cerevisiae.R64-2-1.chrsizes.tsv')
parser.add_argument('-p', dest = 'outpath', type=str, default = 'testout.bed')
args = parser.parse_args()

#import as dataframes, then transform to bedtools object afterwards to verify output
#IT WOULD BE BETTER TO IMPORT AS BEDTOOLS OBJECTS TO VERIFY THE INPUT, BUT I'M LAZY
bed_colnames = ['chrom','start','end','name','score','strand']
bed_coltypes = {'chrom': str, 'start':np.uint64, 'end':np.uint64, 'name': str, 'score':np.float64, 'strand': object}

transcripts = pd.read_table(args.transcripts, names = bed_colnames, dtype = bed_coltypes)
orfs = pd.read_table(args.orfs, names = bed_colnames, dtype = bed_coltypes)
chrom_df = pd.read_table(args.chrsizes, names = ['chrom','size'])
chrom_dict = dict(zip(chrom_df['chrom'].astype(str) ,list(zip(np.zeros(chrom_df.shape[0], dtype=np.int), chrom_df['size']))))

#join transcript and orf info by gene name
indf = transcripts.merge(orfs, how='outer', on='name', suffixes = ['_txn', '_orf'])

#initialize output bed file
outdf = pd.DataFrame(index = np.arange(indf.shape[0]), columns = bed_colnames)

for index, row in indf.iterrows():
    #when a gene is in the transcript and the ORF annotations
    if not np.isnan(row.start_txn) and not np.isnan(row.start_orf):
        if row.strand_txn == "+":
            #check that txn start is left of ORF start
            if row.start_txn <= row.start_orf:
                outdf.iloc[index,] = [row.chrom_txn, max(0, row.start_txn-args.dist), row.start_orf, row['name'], row.score_txn, row.strand_txn]
            else: #if ORF start is left of txn start, FOR NOW WE TAKE THE REGION UPSTREAM OF THE ORF
                outdf.iloc[index,] = [row.chrom_txn, max(0, row.start_orf-args.dist), row.start_orf, row['name'], row.score_txn, row.strand_txn]
        elif row.strand_txn == "-":
            #check that txn start is right of ORF start
            if row.end_orf <= row.end_txn:
                outdf.iloc[index,] = [row.chrom_txn, row.end_orf, row.end_txn+args.dist, row['name'], row.score_txn, row.strand_txn]
            else: #if ORF start is right of txn start, FOR NOW WE TAKE THE REGION UPSTREAM OF THE ORF
                outdf.iloc[index,] = [row.chrom_txn, row.end_orf, row.end_orf+args.dist, row['name'], row.score_txn, row.strand_txn]
    #when a gene is present only in the transcript annotation
    if np.isnan(row.start_orf):
        if row.strand_txn == "+":
            outdf.iloc[index,] = [row.chrom_txn, max(0, row.start_txn-args.dist), row.start_txn+args.dist, row['name'], row.score_txn, row.strand_txn]
        elif row.strand_txn == "-":
            outdf.iloc[index,] = [row.chrom_txn, row.end_txn-args.dist, row.end_txn+args.dist, row['name'], row.score_txn, row.strand_txn]
    #when a gene is present only in the ORF annotation
    if np.isnan(row.start_txn):
        if row.strand_orf == "+":
            outdf.iloc[index,] = [row.chrom_orf, max(0, row.start_orf-args.dist), row.start_orf, row['name'], 0, row.strand_orf]
        elif row.strand_orf == "-":
            outdf.iloc[index,] = [row.chrom_orf, row.end_orf, row.end_orf+args.dist, row['name'], 0, row.strand_orf]

outdf[['start', 'end']] = outdf[['start', 'end']].astype(np.uint64)
outdf['score'] = outdf['score'].astype(np.float64)
outdf = outdf.query('(end-start)>0')
bedout = pybt.BedTool.from_dataframe(outdf).truncate_to_chrom(chrom_dict).remove_invalid().moveto(args.outpath)
