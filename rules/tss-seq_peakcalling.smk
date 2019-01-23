#!/usr/bin/env python

localrules: combine_tss_peaks

rule call_tss_peaks:
    input:
        bw = lambda wc: "coverage/counts/{sample}_tss-seq-counts-SENSE.bw".format(**wc) if wc.species=="experimental" else "coverage/sicounts/{sample}_tss-seq-sicounts-SENSE.bw".format(**wc)
    output:
        smoothed = expand("peakcalling/sample_peaks/{{sample}}_{{species}}-tss-seq-smoothed-bandwidth{bandwidth}-{strand}.bw", strand=["plus","minus"], bandwidth = config["peakcalling"]["bandwidth"]),
        peaks = "peakcalling/sample_peaks/{sample}_{species}-tss-seq-allpeaks.narrowPeak"
    params:
        name = lambda wc: "{sample}_{species}".format(**wc),
        bandwidth = config["peakcalling"]["bandwidth"],
        window = config["peakcalling"]["local_background_window"]
    conda:
        "../envs/peakcalling.yaml"
    log:
        "logs/call_tss_peaks/call_tss_peaks-{sample}-{species}.log"
    shell: """
        (python scripts/tss_peakcalling.py -i {input.bw} -n {params.name} -o peakcalling/sample_peaks -b {params.bandwidth} -w {params.window}) &> {log}
        """

rule tss_peaks_idr:
    input:
        #NOTE: for now we take the first two samples since the IDR script only takes two
        #change this if we find a better way to aggregate results
        lambda wc: ["peakcalling/sample_peaks/" + x + "_{species}-tss-seq-allpeaks.narrowPeak".format(**wc) for x in PASSING if PASSING[x]['group']==wc.group][0:2]
    output:
        allpeaks = "peakcalling/{group}/{group}_{species}-tss-seq-idrpeaks-all.tsv",
        filtered = "peakcalling/{group}/{group}_{species}-tss-seq-idrpeaks-filtered.tsv",
        narrowpeak = "peakcalling/{group}/{group}_{species}-tss-seq-idrpeaks.narrowPeak",
        summits = "peakcalling/{group}/{group}_{species}-tss-seq-idrpeaks-summits.bed",
    params:
        idr = int(-125*log2(config["peakcalling"]["idr"]))
    conda:
        "../envs/peakcalling.yaml"
    log:
        "logs/tss_peaks_idr/tss_peaks_idr-{group}-{species}.log"
    shell: """
        (idr -s {input} --input-file-type narrowPeak --rank q.value -o {output.allpeaks} -l {log} --plot --peak-merge-method max) &> {log}
        (awk '$5>{params.idr} || $9=="inf"' {output.allpeaks} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         tee {output.filtered} | \
         awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $11, $12, $10}}' | \
         sed "s!\(.*\)\(-minus$\|-plus$\)\(.*\)!\1\3!" |
         sort -k1,1 -k2,2n | \
         tee {output.narrowpeak} | \
         awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.summits}) &>> {log}
        """

rule combine_tss_peaks:
    input:
        cond = "peakcalling/{condition}/{condition}_{species}-tss-seq-idrpeaks-filtered.tsv",
        ctrl = "peakcalling/{control}/{control}_{species}-tss-seq-idrpeaks-filtered.tsv",
    output:
        "diff_exp/peaks/{condition}-v-{control}/{condition}-v-{control}_{species}-peaks.bed"
    log:
        "logs/combine_tss_peaks/combine_tss_peaks_{condition}-v-{control}-{species}.log"
    shell: """
        (sort -k1,1 -k2,2n {input.cond} | \
         bedtools multiinter -i stdin <(sort -k1,1 -k2,2n {input.ctrl}) -cluster | \
         cut -f1-3 | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         uniq | \
         awk 'BEGIN{{FS=OFS="\t"}}{{$1 ~ /-plus$/ ? strand="+" : strand="-"; print $1, $2, $3, ".", 0, strand}}' > {output}) &> {log}
        """

