#!/usr/bin/env python

import re
from typing import List
import argparse
import pandas as pd
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy.signal import argrelextrema, gaussian
from scipy.stats import poisson #, entropy
from statsmodels.sandbox.stats.multicomp import multipletests
import pyBigWig as pybw

def gaussian_kernel(sigma: float,
                    width: int,
                    probability: float,
                    position: int,
                    chromsize: int):
    kernel = gaussian(2*width+1,
                      std=sigma)
    #trim to chromosome boundaries
    if position - width < 0:
        kernel = kernel[(width - position):]
    if (position + width) > chromsize:
        if position - width < 0:
            kernel = kernel[:chromsize]
        else:
            kernel = kernel[:(len(kernel) - (position + width - chromsize + 1))]
    #normalize so total density will sum to one
    kernel = np.divide(kernel, kernel.sum())
    kernel = np.multiply(kernel, probability)
    return kernel

def hillclimb(gradient_vector,
              localmax_indices,
              position):
    if position in localmax_indices:
        return position
    elif gradient_vector[position] > 0: #assign nearest local max to the right
        if position+1 in localmax_indices:
            return position+1
        idx = min(np.searchsorted(localmax_indices, position), localmax_indices.size-1)
        return localmax_indices[idx]
    else: #assign nearest local max to the left
        if position-1 in localmax_indices:
            return position-1
        idx = max(np.searchsorted(localmax_indices, position)-1, 0)
        return localmax_indices[idx]

def shorth(array,
           proportion: float):
    sorted_array = np.sort(array)
    width = round(proportion*len(sorted_array))
    if width == len(sorted_array):
        return (sorted_array[0], sorted_array[-1]+1)
    diffs = sorted_array[width:] - sorted_array[:len(sorted_array)-width]
    shorth_start = np.where(diffs == np.min(diffs))[0]
    if len(shorth_start) > 1:
        shorth_start = np.round((np.mean(shorth_start))).astype('uint32')
    shorth_start = np.asscalar(shorth_start)
    return (sorted_array[shorth_start], sorted_array[shorth_start+width]+1)

def trim(series,
         raw_coverage,
         shorth_prop: float):
    data = raw_coverage[series.start:series.end]
    idx = data.nonzero()[0]
    counts = data[idx].astype('uint32')
    positions = np.add(idx, series.start)
    flattened = np.repeat(positions, counts)
    return shorth(flattened, shorth_prop)

def poisson_test(series,
                 raw_coverage,
                 lambda_global: float,
                 lambda_local_window: int):
    logp = 0
    width = series.end-series.start
    if series.chrom.endswith("-plus"):
        window = raw_coverage[max(0, series.start-lambda_local_window) : series.start]
    elif series.chrom.endswith("-minus"):
        window = raw_coverage[series.end : min(len(raw_coverage), series.end+lambda_local_window)]
    if series.counts > 0:
        max_counts = np.max(raw_coverage[series.start : series.end])
        logp += np.log10(poisson.sf(k=max_counts,
                                    mu=max(lambda_global*width, window.sum()/len(window)),
                                    loc=1))
    return np.abs(logp)

# def get_entropy(series,
#                 normalized_coverage):
#     probabilities = normalized_coverage[series.start : series.end]
#     return entropy(probabilities, base=2)

def main(input_paths: List[str],
         output_dir: str,
         sample_names: List[str],
         pilot_bandwidth: int,
         alpha_smooth: float,
         lambda_width: int,
         shorth_prop: float):

    if len(sample_names) != len(input_paths):
        raise ValueError("Number of sample names (-n) does not match \
                number of coverage files provided (-i).")
    if pilot_bandwidth <= 0:
        raise ValueError("Bandwidth for pilot density estimation (-b) \
                must be greater than 0.")
    if not 0 <= alpha_smooth <= 1:
        raise ValueError("Alpha for variable density estimation (-a) must be in \
                [0,1].")
    if lambda_width < 1:
        raise ValueError("Window for finding Poisson lambda (-w) \
                must be positive.")
    if not 0 < shorth_prop <= 1:
        raise ValueError("Proportion for shorth trimming (-f) must be in (0,1].")

    bigwigs = [pybw.open(i) for i in input_paths]

    CHROMS = bigwigs[0].chroms()
    if not all(replicate.chroms() == CHROMS for replicate in bigwigs):
        raise RuntimeError("Not all input bigWigs have the same chromosomes!")
    if not all(chrom.endswith("-plus") or chrom.endswith("-minus") \
            for replicate in bigwigs for chrom in replicate.chroms()):
        raise RuntimeError("""Input bigWigs must be stranded: \
                Coverage for plus and minus strands should be denoted by \
                chromosome names ending in \"-plus\" and \"-minus\".""")

    for replicate_index, bigwig in enumerate(bigwigs):
        df = pd.DataFrame()
        print("processing {sample_name}: {coverage_path}...".format(\
                sample_name=sample_names[replicate_index], \
                coverage_path=input_paths[replicate_index]))

        # open connections to output bigWig files of the smoothed coverage used for peakcalling
        output_bigwigs = {'plus': pybw.open("{directory}/"
                "{sample_name}-tss-seq-smoothed-bandwidth{bandwidth}-plus.bw"
                .format(directory=output_dir, \
                sample_name=sample_names[replicate_index],\
                bandwidth=str(pilot_bandwidth)), "w"),
                          'minus': pybw.open("{directory}/"
                                  "{sample_name}-tss-seq-smoothed-bandwidth{bandwidth}"
                                  "-minus.bw".format(directory=output_dir,\
                                  sample_name=sample_names[replicate_index],\
                                  bandwidth=str(pilot_bandwidth)), "w")}

        output_bigwigs['plus'].addHeader([(re.sub("-plus$", "", chrom), sizes) \
                for chrom, sizes in CHROMS.items() if chrom.endswith("-plus")])
        output_bigwigs['minus'].addHeader([(re.sub("-minus$", "", chrom), sizes) \
                for chrom, sizes in CHROMS.items() if chrom.endswith("-minus")])

        for chrom, size in CHROMS.items():
            raw_coverage = bigwig.values(chrom, 0, size, numpy=True)
            # indices of nonzero datapoints in the coverage
            nonzero_indices = raw_coverage.nonzero()[0]

            # check for empty chromosome
            if not nonzero_indices.size:
                print("no coverage on chromosome: " + chrom)
                continue

            # this can sort of be interpreted as the probability of a start at each base
            normalized_coverage = np.divide(raw_coverage, raw_coverage.sum())

            # stage one: get pilot densities by standard Gaussian KDE
            pilot_density = gaussian_filter1d(normalized_coverage, pilot_bandwidth)

            # stage two: assign adaptive bandwidths
            # g is geometric mean of densities (a scalar)
            g_inv = 1/np.exp(np.dot(normalized_coverage[nonzero_indices],
                                    np.log(pilot_density[nonzero_indices])))
            # the adaptive bandwidth (Gaussian std. dev.) at each nonzero data point
            h_adapt = np.multiply(pilot_bandwidth,
                                  np.power(np.multiply(pilot_density[nonzero_indices],
                                                       g_inv),
                                           -alpha_smooth))

            # get positions within plmin 5 sigma of each nonzero data point
            widths = np.ceil(5*h_adapt).astype('int')
            starts = np.maximum(np.subtract(nonzero_indices, widths), 0)
            ends = np.minimum(np.add(nonzero_indices, widths+1), size)

            # initialize a vector for the smoothed coverage
            smoothed_coverage = np.zeros_like(raw_coverage, dtype='float64')

            # build up smoothed coverage by applying corresponding Gaussian to each data point
            for sigma, width, probability, position, start, end in \
                    zip(h_adapt, widths, normalized_coverage[nonzero_indices], \
                    nonzero_indices, starts, ends):
                smoothed_coverage[start:end] += gaussian_kernel(sigma=sigma,
                                                                width=width,
                                                                probability=probability,
                                                                position=position,
                                                                chromsize=size)

            if chrom.endswith("-plus"):
                output_bigwigs['plus'].addEntries(re.sub("-plus$", "", chrom), \
                        0, values=smoothed_coverage, span=1, step=1)
            elif chrom.endswith("-minus"):
                output_bigwigs['minus'].addEntries(re.sub("-minus$", "", chrom), \
                        0, values=smoothed_coverage, span=1, step=1)

            # assign datapoints to nearest local maximum in the direction of positive derivative
            localmax_indices = argrelextrema(smoothed_coverage, np.greater_equal, order=1)[0]
            localmax_indices = localmax_indices[np.where(smoothed_coverage[localmax_indices] != 0)]

            # deal with flat local maxima (keep the value closest to the middle)
            plateaus = localmax_indices[np.where(np.ediff1d(localmax_indices) == 1)]
            plateaus = np.sort(np.append(plateaus, plateaus+1))
            run = []
            plateaus_grouped = [run]
            expected = None # the expected position of the next local max if in a plateau
            for plateau in plateaus:
                if (plateau == expected) or (expected is None):
                    run.append(plateau)
                else:
                    run = [plateau]
                    plateaus_grouped.append(run)
                expected = plateau + 1
            plateaus_saved = np.array([np.median(group).astype('int') \
                    for group in plateaus_grouped])
            plateaus_dropped = plateaus[np.in1d(plateaus, plateaus_saved, invert=True)]
            localmax_indices = localmax_indices[np.in1d(localmax_indices, \
                    plateaus_dropped, invert=True)]

            gradient = np.gradient(smoothed_coverage)
            hillclimb_vectorized = np.vectorize(hillclimb, \
                    excluded=['gradient_vector', 'localmax_indices'])
            clusters = hillclimb_vectorized(gradient_vector=gradient,
                                            localmax_indices=localmax_indices,
                                            position=nonzero_indices)

            grouped_indices = pd.DataFrame({'cluster':clusters, \
                    'position':nonzero_indices}).groupby('cluster')
            cluster_starts = grouped_indices.apply(lambda g: min(g.position))
            cluster_ends = grouped_indices.apply(lambda g: max(g.position)+1)
            df_chrom = pd.DataFrame({'chrom': chrom,
                                     'start': cluster_starts,
                                     'end': cluster_ends})

            #trim cluster ends to the shorth
            df_chrom[['start', 'end']] = df_chrom.apply(trim,\
                    raw_coverage=raw_coverage,\
                    shorth_prop=shorth_prop,\
                    axis=1).apply(pd.Series, dtype=np.uint32)

            df_chrom['peak'] = df_chrom.apply(lambda x: \
                    np.argmax(raw_coverage[x.start:x.end]), axis=1)
            df_chrom['counts'] = df_chrom.apply(lambda x: \
                    np.sum(raw_coverage[x.start:x.end]), axis=1)

            # calculate Poisson p-value for each cluster
            lambda_chrom = raw_coverage.sum()/size   #chromosome-wide lambda
            df_chrom['logp'] = df_chrom.apply(poisson_test,
                                              raw_coverage=raw_coverage,
                                              lambda_global=lambda_chrom,
                                              lambda_local_window=lambda_width, axis=1)
            # df_chrom['entropy'] = df_chrom.apply(get_entropy,
            #                                      normalized_coverage=normalized_coverage, axis=1)

            df = df.append(df_chrom)
            print(chrom + " complete...")

        #close file connections
        bigwig.close()
        output_bigwigs['plus'].close()
        output_bigwigs['minus'].close()
        df['name'] = "."
        df['strand'] = df["chrom"].apply(lambda x: "+" if x.endswith("-plus") else "-")
        df['logpadj'] = np.abs(np.log10(multipletests(np.power(10, -df['logp']),
                                                      alpha=1e-5,
                                                      method='bonferroni')[1]))
        # df = df[['chrom', 'start', 'end', 'name', 'entropy', \
        df = df[['chrom', 'start', 'end', 'name', 'counts', \
                'strand', 'counts', 'logp', 'logpadj', 'peak']]
        df.to_csv("{directory}/{sample_name}-tss-seq-allpeaks.narrowPeak".\
                format(directory=output_dir, \
                sample_name=sample_names[replicate_index]), \
                sep="\t", header=False, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Call peaks in TSS-seq data.')
    parser.add_argument('-i', dest='input_paths', type=str, nargs='+', \
            help='paths to stranded input coverage files of TSS-seq counts, in bigWig format')
    parser.add_argument('-n', dest='sample_names', type=str, nargs='+', \
            help='sample names, in order of input bigWigs')
    parser.add_argument('-o', dest='output_dir', type=str, default=".", \
            help="directory to write output to")
    parser.add_argument('-b', dest='pilot_bandwidth', type=float, default=10, \
            help='Gaussian kernel bandwidth (i.e., standard deviation) \
            for pilot density estimation')
    parser.add_argument('-a', dest='alpha_smooth', type=float, default=0.2, \
            help='smoothing parameter alpha in [0,1]. \
            Larger values= more flexibility in response to pilot density. \
            Set to 0 to perform standard Gaussian KDE.')
    parser.add_argument('-w', dest='lambda_width', type=int, default=2000, \
            help='width of window in nucleotides for local background Poisson lambda estimation')
    parser.add_argument('-f', dest='shorth_prop', type=float, default=0.95, \
            help='Proportion for shorth trimming, in (0,1]')
    args = parser.parse_args()
    main(input_paths=args.input_paths,
         output_dir=args.output_dir,
         sample_names=args.sample_names,
         pilot_bandwidth=args.pilot_bandwidth,
         alpha_smooth=args.alpha_smooth,
         lambda_width=args.lambda_width,
         shorth_prop=args.shorth_prop)

