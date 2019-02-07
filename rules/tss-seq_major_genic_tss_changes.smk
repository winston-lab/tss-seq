#!/usr/bin/env python

localrules:
    get_gained_utr_bed,
    get_lost_utr_bed,
    get_changed_utr_rna_motifs,
    get_longest_major_utr_sequences,
    find_uorfs,
    get_changed_uorfs

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

rule get_longest_major_utr_sequences:
    input:
        genic_tss_changes = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-nucleotide-changes.tsv",
        orf_anno = os.path.abspath(build_annotations(config["genome"]["orf_annotation"])),
    output:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-longest-5pUTR.bed",
    shell: """
        tail -n +2 {input.genic_tss_changes} | \
        sort -k2,2 | \
        join -1 2 -2 4 -t $'\t' - <(sort -k4,4 {input.orf_anno}) | \
        scripts/get_longest_major_utr.sh | \
        sort -k1,1 -k2,2n > {output}
        """

rule find_uorfs:
    input:
        utrs = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-longest-5pUTR.bed",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"])),
    output:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-uORFs.bed"
    conda:
        "../envs/find_uorfs.yaml"
    shell: """
        python scripts/find_uorfs.py -f {input.fasta} -i {input.utrs} -o {output}
        """

rule get_changed_uorfs:
    input:
        uorfs = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-uORFs.bed",
        utrs = "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-genic-utrs-{change}.bed",
    output:
        "diff_exp/genic-nucleotides/{condition}-v-{control}/{norm}/major-genic-tss-changes/{condition}-v-{control}_tss-seq-{norm}-uORFs-{change}.bed"
    shell: """
        paste {input.uorfs} {input.uorfs} | \
        awk 'BEGIN{{FS=OFS="\t"}} {{$6=="+" ? $3=$2+3 : $2=$3-3; print $0}}' | \
        bedtools intersect -a stdin -b {input.utrs} -s | \
        cut --complement -f1-6 > {output}
        """

