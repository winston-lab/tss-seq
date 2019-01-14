#!/usr/bin/env python

localrules:
    map_counts_to_peaks,
    combine_peak_counts,

rule map_counts_to_peaks:
    input:
        bed = "diff_exp/{condition}-v-{control}/{condition}-v-{control}_{species}-peaks.bed",
        bg = lambda wc: "coverage/counts/{sample}_tss-seq-counts-SENSE.bedgraph".format(**wc) if wc.species=="experimental" else "coverage/sicounts/{sample}_tss-seq-sicounts-SENSE.bedgraph".format(**wc)
    output:
        temp("diff_exp/{condition}-v-{control}/{sample}_tss-seq-{species}-peakcounts.tsv")
    log: "logs/map_counts_to_peaks/map_counts_to_peaks-{condition}-v-{control}_{sample}-{species}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-"$2"-"$3, $4}}' &> {output}) &> {log}
        """

rule combine_peak_counts:
    input:
        lambda wc : ["diff_exp/{condition}-v-{control}/".format(**wc) + x + "_tss-seq-{species}-peakcounts.tsv".format(**wc) for x in get_samples("passing", "libsizenorm", [wc.control, wc.condition])]
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-{species}-peak-counts.tsv"
    params:
        n = lambda wc: 2*len(get_samples("passing", "libsizenorm", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("passing", "libsizenorm", [wc.control, wc.condition]))
    log: "logs/get_peak_counts/get_peak_counts_{condition}-v-{control}_{species}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule differential_expression:
    input:
        expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-experimental-peak-counts.tsv",
        sicounts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-spikein-peak-counts.tsv".format(**wc)
    output:
        results_all = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-all.tsv",
        results_up = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-up.tsv",
        results_down = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-down.tsv",
        results_unch = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-unchanged.tsv",
        normcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-counts-sizefactornorm.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-counts-rlogtransformed.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-qc-plots.svg"
    params:
        samples = lambda wc : get_samples("passing", wc.norm, [wc.control, wc.condition]),
        groups = lambda wc : [PASSING[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition])],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    conda: "../envs/diff_exp.yaml"
    script:
        "../scripts/call_diffexp_peaks.R"

rule diffexp_results_to_narrowpeak:
    input:
        condition_coverage = lambda wc: expand("coverage/{norm}/{sample}_tss-seq-{norm}-SENSE.bw", sample=get_samples("passing", wc.norm, wc.condition), norm=wc.norm),
        control_coverage = lambda wc: expand("coverage/{norm}/{sample}_tss-seq-{norm}-SENSE.bw", sample=get_samples("passing", wc.norm, wc.control), norm=wc.norm),
        diffexp_results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        summit_bed = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}-summits.bed",
    conda: "../envs/peakcalling.yaml"
    log: "logs/diffexp_results_to_narrowpeak/diffexp_results_to_narrowpeak-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/diffexp_results_to_narrowpeak.py -i {input.condition_coverage} -j {input.control_coverage} -d {input.diffexp_results} -n {output.narrowpeak} -b {output.summit_bed}) &> {log}
        """

rule summarise_diffexp_results:
    input:
        total = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-all.tsv",
        genic = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-all.tsv",
        intragenic = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-all.tsv",
        antisense = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-antisense-all.tsv",
        convergent = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-all.tsv",
        divergent = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-all.tsv",
        intergenic = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intergenic-all.tsv",
    output:
        summary_table = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-summary.tsv",
        mosaic = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-mosaic.svg",
        maplot = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-maplot.svg",
        volcano = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-volcano.svg",
        volcano_free = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-volcano-freescale.svg",
    params:
        lfc = config["deseq"]["fold-change-threshold"],
        alpha = config["deseq"]["fdr"]
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_diffexp_summary.R"

