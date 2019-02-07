#!/usr/bin/env python

localrules:
    get_gained_utr_bed,
    get_lost_utr_bed,
    get_changed_utr_rna_motifs

rule get_major_genic_tss_changes:
    input:
        results = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-genic-nucleotides-diffexp-results-all.tsv",
    output:
        tsv = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-nucleotide-changes.tsv",
        expression_plot = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-nucleotide-expression-plot.svg",
        shift_plot = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-nucleotide-shift-plot.svg",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/genic_nucleotide_changes.R"

rule get_gained_utr_bed:
    input:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-nucleotide-changes.tsv",
    output:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-utrs-gained.bed",
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} NR>1 && $6<-6 {{if ($3=="-") {{start=$4+1; end=$5+1}} else {{start=$5; end=$4}}; print $1, start, end, $2, $6, $3}}' {input} > {output}
        """

rule get_lost_utr_bed:
    input:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-nucleotide-changes.tsv",
    output:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-utrs-lost.bed",
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}} NR>1 && $6>6 {{if ($3=="-") {{end=$4+1; start=$5+1}} else {{end=$5; start=$4}}; print $1, start, end, $2, $6, $3}}' {input} > {output}
        """

rule get_changed_utr_rna_motifs:
    input:
        utrs = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-utrs-{change}.bed",
        motifs = build_annotations("motifs/" + config["genome"]["name"] + "_all_rna_motifs.bed")
    output:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-utrs-{change}-rna-motifs.tsv"
    shell: """
        bedtools intersect -a {input.utrs} -b {input.motifs} -s -F 1 > {output}
        """






