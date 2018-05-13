#!/usr/bin/env python

rule map_counts_to_peaks:
    input:
        bed = "diff_exp/{condition}-v-{control}/{condition}-v-{control}_{type}-peaks.bed",
        bg = lambda wc: "coverage/counts/" + wc.sample + "_tss-seq-counts-SENSE.bedgraph" if wc.type=="experimental" else "coverage/sicounts/" + wc.sample + "-tss-sicounts-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/{sample}_tss-seq-{type}-peakcounts.tsv")
    log: "logs/map_counts_to_peaks/map_counts_to_peaks-{condition}-v-{control}_{sample}-{type}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-"$2"-"$3, $4}}' &> {output}) &> {log}
        """

rule combine_peak_counts:
    input:
        lambda wc : ["diff_exp/" + wc.condition + "-v-" + wc.control + "/" + x + "_tss-seq-" + wc.type + "-peakcounts.tsv" for x in getsamples(wc.control, wc.condition)]
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-{type}-peak-counts.tsv"
    params:
        n = lambda wc: 2*len(getsamples(wc.control, wc.condition)),
        names = lambda wc: "\t".join(getsamples(wc.control, wc.condition))
    log: "logs/get_peak_counts/get_peak_counts_{condition}-v-{control}_{type}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule differential_expression:
    input:
        expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-experimental-peak-counts.tsv",
        sicounts = lambda wc: "diff_exp/" + wc.condition + "-v-" + wc.control + "/" + wc.condition + "-v-" + wc.control + "_tss-seq-spikein-peak-counts.tsv" if wc.norm=="spikenorm" else []
    output:
        results_all = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-all.tsv",
        results_up = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-up.tsv",
        results_down = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-down.tsv",
        results_unch = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-unchanged.tsv",
        normcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-counts-sizefactornorm.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-counts-rlogtransformed.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-qc-plots.svg"
    params:
        samples = lambda wc : getsamples(wc.control, wc.condition),
        groups = lambda wc : [PASSING[x]["group"] for x in getsamples(wc.control, wc.condition)],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    script:
        "scripts/call_de_peaks.R"

rule diffexp_results_to_narrowpeak:
    input:
        condition_coverage = lambda wc: expand("coverage/{norm}/{sample}_tss-seq-{norm}-SENSE.bw", sample=getsamples(wc.condition, wc.condition), norm=wc.norm),
        control_coverage = lambda wc: expand("coverage/{norm}/{sample}_tss-seq-{norm}-SENSE.bw", sample=getsamples(wc.control, wc.control), norm=wc.norm),
        diffexp_results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        summit_bed = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}-summits.bed",
    log: "logs/diffexp_results_to_narrowpeak/diffexp_results_to_narrowpeak-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/diffexp_results_to_narrowpeak.py -i {input.condition_coverage} -j {input.control_coverage} -d {input.diffexp_results} -n {output.narrowpeak} -b {output.summit_bed}) &> {log}
        """

#TODO: account for double-counted peaks when a peak overlaps more than one annotation (more than one genic region, for example)
rule summarise_de_results:
    input:
        total = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-all.tsv",
        genic = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}-results-{norm}-all-genic.tsv",
        intragenic = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-all-intragenic.tsv",
        antisense = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}-results-{norm}-all-antisense.tsv",
        convergent = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}-results-{norm}-all-convergent.tsv",
        divergent = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}-results-{norm}-all-divergent.tsv",
        intergenic = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}-results-{norm}-all-intergenic.tsv",
    params:
        lfc = config["deseq"]["fold-change-threshold"],
        alpha = config["deseq"]["fdr"]
    output:
        summary_table = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-summary.tsv",
        summary = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-summary.svg",
        maplot = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-maplot.svg",
        volcano = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-volcano.svg",
        volcano_free = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-volcano-freescale.svg",
    script: "scripts/de_summary.R"

