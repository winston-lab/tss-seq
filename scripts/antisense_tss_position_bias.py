#!/usr/bin/env python

import argparse
import pandas as pd
import pybedtools as pybt
from collections import Counter

parser = argparse.ArgumentParser(description='Define possible positions for antisense TSSs to occur.')
parser.add_argument('-d', dest = 'diffexp_path', type=str)
parser.add_argument('-c', dest = 'chrom_sizes', type=str, nargs="+")
parser.add_argument('-a', dest = 'tss_path', type=str)
parser.add_argument('-r', dest = 'rel_path', type=str)
parser.add_argument('-s', dest = 'cps_path', type=str)
args = parser.parse_args()

def main():
    diffexp_results = pd.read_csv(args.diffexp_path, sep="\t")
    chrom_sizes = {args.chrom_sizes[i]:(0, int(args.chrom_sizes[i+1])) for i in range(0, len(args.chrom_sizes), 2)}

    tss_values = Counter()
    cps_values = Counter()
    relative_values = Counter()

    # for each transcript+peak pair:
        # extend the transcript by (peak length-1) on both sides and switch strand to find all possible antisense positions
        # shift annotation by the summit to get possible summit positions

    for index, row in diffexp_results.iterrows():
        peak_length = row['end']-row['start']

        bed = pybt.BedTool("\t".join([str(row['transcript_chrom']), str(row['transcript_start']), str(row['transcript_end']), str(row['transcript_name']), str(row['transcript_score']), str(row['strand'])]), from_string=True).slop(b=(peak_length-1), g=chrom_sizes).shift(p=row['peak_summit'], m=-row['peak_summit'], g=chrom_sizes)

        for interval in bed:
            if interval.strand=="+":
                tss_values.update(range(interval.start-row['transcript_start'], interval.end-row['transcript_start']))
                cps_values.update(range(interval.start-row['transcript_end'], interval.end-row['transcript_end']))
                relative_values.update([round(x/(row['transcript_end']-row['transcript_start']), 3) for x in range(interval.start-row['transcript_start'], interval.end-row['transcript_start'])])
            else:
                tss_values.update(range(row['transcript_end']-interval.start, row['transcript_end']-interval.end))
                cps_values.update(range(row['transcript_start']-interval.start, row['transcript_start']-interval.end))
                relative_values.update([round(x/(row['transcript_end']-row['transcript_start']), 3) for x in range(row['transcript_end']-interval.start, row['transcript_end']-interval.end)])


    tss_df = pd.DataFrame.from_dict(tss_values, orient="index", columns=["count"]).reset_index().sort_values(["index"])
    cps_df = pd.DataFrame.from_dict(cps_values, orient="index", columns=["count"]).reset_index().sort_values(["index"])
    relative_df = pd.DataFrame.from_dict(relative_values, orient="index", columns=["count"]).reset_index().sort_values(["index"])

    tss_df.to_csv(args.tss_path, sep="\t", index=False)
    cps_df.to_csv(args.cps_path, sep="\t", index=False)
    relative_df.to_csv(args.rel_path, sep="\t", index=False)

if __name__ == '__main__':
    main()

