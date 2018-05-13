#!/usr/bin/env python

rule call_tss_peaks:
    input:
        bw = lambda wc: "coverage/counts/" + wc.sample + "_tss-seq-counts-SENSE.bw" if wc.type=="experimental" else "coverage/sicounts/" + wc.sample + "_tss-seq-sicounts-SENSE.bw"
    output:
        smoothed = expand("peakcalling/sample_peaks/{{sample}}_{{type}}-smoothed-bw{bandwidth}-{strand}.bw", strand=["plus","minus"], bandwidth = config["peakcalling"]["bandwidth"]),
        peaks = "peakcalling/sample_peaks/{sample}_{type}-allpeaks.narrowPeak"
    params:
        name = lambda wc: wc.sample + "_" + wc.type,
        bandwidth = config["peakcalling"]["bandwidth"],
        window = config["peakcalling"]["local-bg-window"]
    log: "logs/call_tss_peaks/call_tss_peaks-{sample}-{type}.log"
    shell: """
        (python scripts/tss-peakcalling.py -i {input.bw} -n {params.name} -w {params.window} -b {params.bandwidth} -o peakcalling/sample_peaks) &> {log}
        """

rule tss_peaks_idr:
    input:
        #NOTE: for now we take the first two samples since the IDR script only takes two
        #change this if we find a better way to aggregate results
        lambda wc: ["peakcalling/sample_peaks/" + x + "_" + wc.type + "-allpeaks.narrowPeak" for x in PASSING if PASSING[x]['group']==wc.group][0:2]
    output:
        allpeaks = "peakcalling/{group}/{group}_{type}-idrpeaks-all.tsv",
        filtered = "peakcalling/{group}/{group}_{type}-idrpeaks-filtered.tsv",
        narrowpeak = "peakcalling/{group}/{group}_{type}-idrpeaks.narrowPeak",
        summits = "peakcalling/{group}/{group}_{type}-idrpeaks-summits.bed",
    params:
        idr = int(-125*log2(config["peakcalling"]["idr"]))
    log: "logs/tss_peaks_idr/tss_peaks_idr-{group}-{type}.log"
    shell: """
        idr -s {input} --input-file-type narrowPeak --rank q.value -o {output.allpeaks} -l {log} --plot --peak-merge-method max
        (awk '$5>{params.idr} || $9=="inf"' {output.allpeaks} | tee {output.filtered} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $11, $12, $10}}' | sed "s/-minus//g;s/-plus//g" | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.summits}) &> {log}
        """

rule combine_tss_peaks:
    input:
        cond = "peakcalling/{condition}/{condition}_{type}-idrpeaks-filtered.tsv",
        ctrl = "peakcalling/{control}/{control}_{type}-idrpeaks-filtered.tsv",
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}_{type}-peaks.bed"
    shell: """
        sort -k1,1 -k2,2n {input.cond} | bedtools multiinter -i stdin <(sort -k1,1 -k2,2n {input.ctrl}) -cluster | cut -f1-3 | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

