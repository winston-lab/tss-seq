#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gsmooth
from scipy.signal import argrelmax
from scipy.stats import poisson
from statsmodels.sandbox.stats.multicomp import multipletests
import pyBigWig as pybw

parser = argparse.ArgumentParser(description='Smooth bigwig file with Gaussian kernel of given bandwidth.')
parser.add_argument('-i', dest = 'infiles', type=str, nargs='+', help='input BigWig files')
parser.add_argument('-b', dest = 'bandwidth', type=int, default = 20, help='Gaussian kernel bandwidth (standard deviation)')
parser.add_argument('-o', dest = 'out', type=str, default=".", help="output files destination")
parser.add_argument('-n', dest = 'names', type=str, nargs='+', help='sample names, in order of input bigwigs')
parser.add_argument('-w', dest = 'wwidth', type=int, default=2000, help='half-width of window for local background')
args = parser.parse_args()

def poistest(series, rr, ll, ws):
    logp = 0
    width = series.end-series.start

    window = rr[max(0, series.start-ws):min(len(rr), series.end+ws)]
    counts = rr[series.start:series.end].sum()

    if counts>0:
        logp += np.log10(poisson.sf(k=counts, mu=max(ll*width, window.sum()/len(window)), loc=1))

    return np.abs(logp)

#open connections to input bigwigs
inbw = [pybw.open(i) for i in args.infiles]
nreps = len(args.infiles)

#verify bigwigs are from same genome
CHROMS = inbw[0].chroms()
if not all(replicate.chroms() == CHROMS for replicate in inbw):
    raise RuntimeError('Not all input bigwigs are from the same genome!')

for rep, bw in enumerate(inbw):

    df = pd.DataFrame()
    print("processing " + args.names[rep] + ": " + args.infiles[rep] + "..." )

    #open connections to output bigwigs
    outbw = {'plus': pybw.open(args.out + "/" + args.names[rep] + "-smoothed-bw" + str(args.bandwidth) + "-plus.bw", "w"),
             'minus':pybw.open(args.out + "/" + args.names[rep] + "-smoothed-bw" + str(args.bandwidth) + "-minus.bw", "w")}
    outbw['plus'].addHeader([(c.replace("-plus", ""), ss) for c, ss in CHROMS.items() if "-plus" in c])
    outbw['minus'].addHeader([(c.replace("-minus", ""), ss) for c, ss in CHROMS.items() if "-minus" in c])

    for chrom, size in CHROMS.items():
        #import data
        raw = bw.values(chrom, 0, size, numpy=True)

        #smooth by Gaussian convolution
        smoothed = gsmooth(raw, sigma=args.bandwidth, order=0, mode='mirror')
        if "plus" in chrom:
            outbw['plus'].addEntries(chrom.replace("-plus", ""), 0, values=smoothed, span=1, step=1)
        elif "minus" in chrom:
            outbw['minus'].addEntries(chrom.replace("-minus", ""), 0, values=smoothed, span=1, step=1)

        #find location of local maxima
        locmaxidx = argrelmax(smoothed)[0]

        #get nonzero datapoints
        nzidx = raw.nonzero()[0]

        #filter out nonzero datapoints that are not above background
        #background is Poisson with lambda in local window or across chromosome (a la MACS)
        wsize = args.wwidth
        lbdchrom = raw.sum()/size   #chromosome-wide lambda
        lbdlocal = np.array([max(lbdchrom, raw[max(0,pos-wsize):min(size,pos+wsize)].sum() / (min(size,pos+wsize)-max(0,pos-wsize))) for pos in nzidx])
        bgprob = poisson.sf(raw[nzidx], lbdlocal, loc=1)
        nzidxfiltered  = nzidx[bgprob<0.01]

        #assign datapoints to nearest local maximum, take min and max positions as cluster boundaries
        #NOTE: is this right or do we have to iterate as in mean-shift? seems too simple but then again it is 1D...
        nearestlocmax = [locmaxidx[np.abs(np.subtract(locmaxidx,pos)).argmin()] for pos in nzidxfiltered]
        grouped = pd.DataFrame({'max':nearestlocmax, 'pos':nzidxfiltered}).groupby('max')

        starts = grouped.apply(lambda g: min(g.pos))
        #check for empty dataframe
        if starts.empty:
            print(chrom)
            break
        ends = grouped.apply(lambda g: max(g.pos)+1)

        chromdf = pd.DataFrame({'chrom':chrom, 'start':starts, 'end':ends})
        chromdf['peak'] = chromdf.apply(lambda x: np.argmax(raw[x.start:x.end]), axis=1)
        chromdf['counts'] = chromdf.apply(lambda x: np.sum(raw[x.start:x.end]), axis=1)

        #calculate Poisson p-value for each cluster
        chromdf['logp'] = chromdf.apply(poistest, rr = raw, ll = lbdchrom, ws=wsize, axis=1)

        df = df.append(chromdf)
        print(chrom)

    #close file connections
    bw.close()
    outbw['plus'].close()
    outbw['minus'].close()

    df['name'] = "."
    df['strand'] = df["chrom"].apply(lambda x: "+" if "plus" in x else "-")
    df['logpadj'] = np.abs(np.log10(multipletests(np.power(10, -df['logp']), alpha=1e-5, method='bonferroni')[1]))
    df = df[['chrom', 'start', 'end', 'name', 'counts', 'strand', 'counts','logp', 'logpadj', 'peak']]
    df.to_csv(args.out + "/" + args.names[rep] + "-allpeaks.narrowPeak" , sep="\t", header=False, index=False)
