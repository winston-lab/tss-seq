#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
# from scipy.signal import argrelmax, gaussian
from scipy.signal import argrelextrema, gaussian
from scipy.stats import poisson
from statsmodels.sandbox.stats.multicomp import multipletests
import pyBigWig as pybw

parser = argparse.ArgumentParser(description='Smooth bigwig file with Gaussian kernel of given bandwidth.')
parser.add_argument('-i', dest = 'infiles', type=str, nargs='+', help='input BigWig files')
parser.add_argument('-b', dest = 'bandwidth', type=int, default = 10, help='Gaussian kernel bandwidth (standard deviation) for pilot density estimation')
parser.add_argument('-o', dest = 'out', type=str, default=".", help="output files destination")
parser.add_argument('-n', dest = 'names', type=str, nargs='+', help='sample names, in order of input bigwigs')
parser.add_argument('-w', dest = 'wwidth', type=int, default=2000, help='width of window for local background')
parser.add_argument('-a', dest = 'alpha', type=float, default=0.2, help='smoothing parameter alpha in (0,1). Larger values= more flexibility in response to pilot density.')
parser.add_argument('-f', dest = 'frac', type=float, default=0.95, help='fraction for shorth trimming')
args = parser.parse_args()

def gkernel(s, width, prob, pos, size):
    out = gaussian(2*width+1, std=s)
    #trim to chromosome boundaries
    if pos-width<0:
        out = out[width-pos:]
    if (pos+width)>size:
        if pos-width<0:
            out = out[:size]
        else:
            out = out[:len(out)-(width+pos-size+1)]
    #normalize so total density will sum to one
    out = np.divide(out, out.sum())
    out = np.multiply(out, prob)
    return out

def hillclimb(grad, locmaxidx, pos):
    if pos in locmaxidx:
        return pos
    elif grad[pos] > 0: #assign nearest local max to the right
        if pos+1 in locmaxidx:
            return pos+1
        idx = min(np.searchsorted(locmaxidx, pos), locmaxidx.size-1)
        return locmaxidx[idx]
    else: #assign nearest local max to the left
        if pos-1 in locmaxidx:
            return pos-1
        idx = max(np.searchsorted(locmaxidx, pos)-1, 0)
        return locmaxidx[idx]

def poistest(series, rr, ll, ws):
    logp = 0
    width = series.end-series.start
    if 'plus' in series.chrom:
        window = rr[max(0, series.start-ws):series.start]
    elif 'minus' in series.chrom:
        window = rr[series.end:min(len(rr), series.end+ws)]
    if series.counts>0:
        logp += np.log10(poisson.sf(k=series.counts, mu=max(ll*width, window.sum()/len(window)), loc=1))
    return np.abs(logp)

def shorth(array, frac):
    sx = np.sort(array)
    width = round(frac*len(sx))
    if width==len(sx):
        return (sx[0], sx[-1]+1)
    diffs = sx[width:] - sx[:len(sx)-width]
    q = np.where(diffs==np.min(diffs))[0]
    if len(q)>1:
        q = np.round((np.mean(q))).astype('uint32')
    q = np.asscalar(q)
    return (sx[q],sx[q+width]+1)

def trim(series, rr, frac):
    data = rr[series.start:series.end]
    idx = data.nonzero()[0]
    counts = data[idx].astype('uint32')
    positions = np.add(idx, series.start)
    onedim = np.repeat(positions, counts)
    return shorth(onedim, frac)

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
        nzidx = raw.nonzero()[0] #indices of nonzero datapoints

        #check for empty chromosome
        if not nzidx.size:
            print("no coverage on chromosome: " + chrom)
            continue

        normalized = np.divide(raw, raw.sum()) #this can be interpreted as the probability at each base

        #stage one: get pilot densities
        pilot = gaussian_filter1d(normalized, args.bandwidth)

        #stage two: assign adaptive bandwidths
        ginv = 1/np.exp(np.dot(normalized[nzidx], np.log(pilot[nzidx]))) #g is geometric mean of densities
        h_adapt = np.multiply(args.bandwidth, np.power(np.multiply(pilot[nzidx], ginv), -args.alpha))

        w = np.ceil(5*h_adapt).astype('int') #discretized 5 sigma cutoff
        starts = np.maximum(np.subtract(nzidx, w), 0)
        ends = np.minimum(np.add(nzidx, w+1), size)

        smoothed = np.zeros_like(raw, dtype='float64')

        for sigma, width, prob, pos, start, end in zip(h_adapt, w, normalized[nzidx], nzidx, starts, ends):
            smoothed[start:end] += gkernel(sigma, width, prob, pos, size)

        if "plus" in chrom:
            outbw['plus'].addEntries(chrom.replace("-plus", ""), 0, values=smoothed, span=1, step=1)
        elif "minus" in chrom:
            outbw['minus'].addEntries(chrom.replace("-minus", ""), 0, values=smoothed, span=1, step=1)

        # assign datapoints to nearest local maximum in the direction of positive derivative
        locmaxidx = argrelextrema(smoothed, np.greater_equal, order=1)[0]
        locmaxidx = locmaxidx[np.where(smoothed[locmaxidx] != 0)]

        #deal with flat local maxima (keep the value closest to the middle)
        flatmaxima = locmaxidx[np.where(np.ediff1d(locmaxidx)==1)]
        flatmaxima = np.sort(np.append(flatmaxima, flatmaxima+1))
        run = []
        groupedflatmaxima = [run]
        expect=None
        for v in flatmaxima:
            if (v==expect) or (expect is None):
                run.append(v)
            else:
                run = [v]
                groupedflatmaxima.append(run)
            expect = v+1
        keepflatmaxima = np.array([np.median(i).astype('int') for i in groupedflatmaxima])
        dropflatmaxima = flatmaxima[np.in1d(flatmaxima, keepflatmaxima, invert=True)]
        locmaxidx = locmaxidx[np.in1d(locmaxidx, dropflatmaxima, invert=True)]

        grad = np.gradient(smoothed)
        vhillclimb = np.vectorize(hillclimb, excluded=['grad', 'locmaxidx'])
        clusters = vhillclimb(grad=grad, locmaxidx=locmaxidx, pos=nzidx)

        grouped = pd.DataFrame({'cluster':clusters, 'position':nzidx}).groupby('cluster')
        clust_starts = grouped.apply(lambda g: min(g.position))
        clust_ends = grouped.apply(lambda g: max(g.position)+1)
        chromdf = pd.DataFrame({'chrom':chrom, 'start':clust_starts, 'end':clust_ends})

        #trim cluster ends to the shorth
        chromdf[['start','end']] = chromdf.apply(trim, rr=raw, frac=args.frac, axis=1).apply(pd.Series, dtype=np.uint32)

        chromdf['peak'] = chromdf.apply(lambda x: np.argmax(raw[x.start:x.end]), axis=1)
        chromdf['counts'] = chromdf.apply(lambda x: np.sum(raw[x.start:x.end]), axis=1)

        # calculate Poisson p-value for each cluster
        lbdchrom = raw.sum()/size   #chromosome-wide lambda
        chromdf['logp'] = chromdf.apply(poistest, rr = raw, ll = lbdchrom, ws=args.wwidth, axis=1)

        df = df.append(chromdf)
        print(chrom + " complete...")

    #close file connections
    bw.close()
    outbw['plus'].close()
    outbw['minus'].close()

df['name'] = "."
df['strand'] = df["chrom"].apply(lambda x: "+" if "plus" in x else "-")
df['logpadj'] = np.abs(np.log10(multipletests(np.power(10, -df['logp']), alpha=1e-5, method='bonferroni')[1]))
df = df[['chrom', 'start', 'end', 'name', 'counts', 'strand', 'counts','logp', 'logpadj', 'peak']]
df.to_csv(args.out + "/" + args.names[rep] + "-allpeaks.narrowPeak" , sep="\t", header=False, index=False)
