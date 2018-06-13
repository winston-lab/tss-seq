#!/usr/bin/env python

import argparse
import pandas as pd
import pybedtools as pybt
from collections import Counter

parser = argparse.ArgumentParser(description='Define possible positions for intragenic TSSs to occur.')
parser.add_argument('-d', dest = 'diffexp_path', type=str)
parser.add_argument('-c', dest = 'chrom_sizes', type=str)
parser.add_argument('-g', dest = 'genic_region_path', type=str)
parser.add_argument('-a', dest = 'atg_path', type=str)
parser.add_argument('-r', dest = 'rel_path', type=str)
parser.add_argument('-s', dest = 'stop_path', type=str)
args = parser.parse_args()

def main():
    diffexp_results = pd.read_csv(args.diffexp_path, sep="\t")
    genic_regions = pybt.BedTool(args.genic_region_path)

    atg_values = Counter()
    stop_values = Counter()
    relative_values = Counter()

    # for each ORF+peak pair:
        # extend the ORF by (peak length-1) on both sides to find all possible intragenic positions
        # subtract regions which overlap genic regions
        # shift annotation by the summit to get possible summit positions

    for index, row in diffexp_results.iterrows():
        peak_length = row['end']-row['start']

        bed = pybt.BedTool("\t".join([row['orf_chrom'], str(row['orf_start']), str(row['orf_end']), row['orf_name'], str(row['orf_score']), row['orf_strand']]), from_string=True).slop(b=(peak_length-1), g=args.chrom_sizes).subtract(genic_regions, s=True).shift(p=row['peak_summit'], m=-row['peak_summit'], g=args.chrom_sizes)

        for interval in bed:
            if interval.strand=="+":
                atg_values.update(range(interval.start-row['orf_start'], interval.end-row['orf_start']))
                stop_values.update(range(interval.start-row['orf_end'], interval.end-row['orf_end']))
                relative_values.update([round(x/(row['orf_end']-row['orf_start']), 3) for x in range(interval.start-row['orf_start'], interval.end-row['orf_start'])])
            else:
                atg_values.update(range(row['orf_end']-interval.start, row['orf_end']-interval.end))
                stop_values.update(range(row['orf_start']-interval.start, row['orf_start']-interval.end))
                relative_values.update([round(x/(row['orf_end']-row['orf_start']), 3) for x in range(row['orf_end']-interval.start, row['orf_end']-interval.end)])


    atg_df = pd.DataFrame.from_dict(atg_values, orient="index", columns=["count"]).reset_index().sort_values(["index"])
    stop_df = pd.DataFrame.from_dict(stop_values, orient="index", columns=["count"]).reset_index().sort_values(["index"])
    relative_df = pd.DataFrame.from_dict(relative_values, orient="index", columns=["count"]).reset_index().sort_values(["index"])

    atg_df.to_csv(args.atg_path, sep="\t", index=False)
    stop_df.to_csv(args.stop_path, sep="\t", index=False)
    relative_df.to_csv(args.rel_path, sep="\t", index=False)

if __name__ == '__main__':
    main()

