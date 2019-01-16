#!/usr/bin/env python

localrules:
    map_counts_to_annotation,
    combine_annotation_counts,
    aggregate_genic_tss

rule aggregate_genic_tss:
    input:
        "diff_exp/peaks/{condition}-v-{control}/libsizenorm/genic/{condition}-v-{control}_tss-seq-libsizenorm-peaks-diffexp-results-genic-all.tsv",
    output:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{condition}-v-{control}_{species}-genic-nucleotides.bed"
    shell: """
        tail -n +2 {input} | \
        cut -f1-6,19 | \
        sort -k7,7 | \
        bedtools groupby -g 7 -c 1,2,3,6 -o first,min,max,first | \
        awk 'BEGIN{{FS=OFS="\t"}}{{$5=="+" ? strand="-plus" : strand="-minus"; print $2strand, $3, $4, $1, 0, $5}}' | \
        bedtools makewindows -b stdin -w 1 -i src | \
        LC_COLLATE=C sort -k1,1 -k2,2n | \
        awk 'BEGIN{{FS=OFS="\t"}}{{$1 ~ /-plus$/ ? strand="+" : strand="-"; print $1, $2, $3, $4, 0, strand}}' > {output}
        """

rule map_counts_to_annotation:
    input:
        bed = "diff_exp/{annotation}/{condition}-v-{control}/{condition}-v-{control}_{species}-{annotation}.bed",
        bg = lambda wc: "coverage/counts/{sample}_tss-seq-counts-SENSE.bedgraph".format(**wc) if wc.species=="experimental" else "coverage/sicounts/{sample}_tss-seq-sicounts-SENSE.bedgraph".format(**wc)
    output:
        temp("diff_exp/{annotation}/{condition}-v-{control}/{sample}_tss-seq-{species}-counts-{annotation}.tsv")
    log:
        "logs/map_counts_to_annotation/map_counts_to_annotation-{condition}-v-{control}_{sample}-{species}-{annotation}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule combine_annotation_counts:
    input:
        lambda wc : ["diff_exp/{annotation}/{condition}-v-{control}/".format(**wc) + x + "_tss-seq-{species}-counts-{annotation}.tsv".format(**wc) for x in get_samples("passing", "libsizenorm", [wc.control, wc.condition])]
    output:
        "diff_exp/{annotation}/{condition}-v-{control}/{condition}-v-{control}_tss-seq-{species}-counts-{annotation}.tsv"
    params:
        n = lambda wc: 7*len(get_samples("passing", "libsizenorm", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("passing", "libsizenorm", [wc.control, wc.condition]))
    log:
        "logs/get_peak_counts/get_peak_counts_{condition}-v-{control}_{species}-{annotation}.log"
    shell: """
        (paste {input} | \
         cut -f$(paste -d, <(echo "1-6") <(seq -s, 7 7 {params.n})) | \
         cat <(echo -e "chrom\tstart\tend\tname\tscore\tstrand\t{params.names}" ) - > {output}) &> {log}
        """

rule differential_expression:
    input:
        exp_counts = "diff_exp/{annotation}/{condition}-v-{control}/{condition}-v-{control}_tss-seq-experimental-counts-{annotation}.tsv",
        spike_counts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_exp/peaks/{condition}-v-{control}/{condition}-v-{control}_tss-seq-spikein-counts-peaks.tsv".format(**wc)
    output:
        results_all = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{annotation}-diffexp-results-all.tsv",
        results_up = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{annotation}-diffexp-results-up.tsv",
        results_down = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{annotation}-diffexp-results-down.tsv",
        results_unchanged = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{annotation}-diffexp-results-unchanged.tsv",
        counts_norm = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{annotation}-counts-sizefactornorm.tsv",
        counts_rlog = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{annotation}-counts-rlogtransformed.tsv",
        qc_plots = "diff_exp/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{annotation}-qc-plots.svg"
    params:
        samples = lambda wc : get_samples("passing", wc.norm, [wc.control, wc.condition]),
        groups = lambda wc : [PASSING[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition])],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    conda:
        "../envs/diff_exp.yaml"
    script:
        "../scripts/differential_expression.R"


rule diffexp_results_to_narrowpeak:
    input:
        condition_coverage = lambda wc: expand("coverage/{norm}/{sample}_tss-seq-{norm}-SENSE.bw", sample=get_samples("passing", wc.norm, wc.condition), norm=wc.norm),
        control_coverage = lambda wc: expand("coverage/{norm}/{sample}_tss-seq-{norm}-SENSE.bw", sample=get_samples("passing", wc.norm, wc.control), norm=wc.norm),
        diffexp_results = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-{direction}.tsv",
    output:
        narrowpeak = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-{direction}.narrowpeak",
        summit_bed = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-{direction}-summits.bed",
    conda:
        "../envs/peakcalling.yaml"
    log:
        "logs/diffexp_results_to_narrowpeak/diffexp_results_to_narrowpeak-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/diffexp_results_to_narrowpeak.py -i {input.condition_coverage} -j {input.control_coverage} -d {input.diffexp_results} -n {output.narrowpeak} -b {output.summit_bed}) &> {log}
        """

rule summarise_diffexp_results:
    input:
        total = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-all.tsv",
        genic = "diff_exp/peaks/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-genic-all.tsv",
        intragenic = "diff_exp/peaks/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-intragenic-all.tsv",
        antisense = "diff_exp/peaks/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-antisense-all.tsv",
        convergent = "diff_exp/peaks/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-convergent-all.tsv",
        divergent = "diff_exp/peaks/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-divergent-all.tsv",
        intergenic = "diff_exp/peaks/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-intergenic-all.tsv",
    output:
        summary_table = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-summary.tsv",
        mosaic = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-mosaic.svg",
        maplot = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-maplot.svg",
        volcano = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-volcano.svg",
        volcano_free = "diff_exp/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-volcano-freescale.svg",
    params:
        lfc = config["deseq"]["fold-change-threshold"],
        alpha = config["deseq"]["fdr"]
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_diffexp_summary.R"

