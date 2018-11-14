#!/usr/bin/env python

localrules:
    get_random_motifs,
    get_meme_sequences

#bedtools intersect peaks with fimo motifs
#0. with the summit of the peak as reference, extend annotation to upstream and 'downstream' distances
#1. merge overlapping (but not book-ended) features
#2. intersect with motif file
rule get_overlapping_motifs:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.narrowpeak",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        motifs = build_annotations("motifs/" + config["genome"]["name"] + "_allmotifs.bed")
    output:
        "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-allFIMOresults.tsv.gz",
    params:
        upstr = config["motifs"]["enrichment-upstream"],
        dnstr = config["motifs"]["enrichment-downstream"]
    log: "logs/get_upstream_motifs/get_upstream_motifs-{condition}-v-{control}-{norm}-{category}-{direction}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.peaks} | \
        bedtools slop -l {params.upstr} -r {params.dnstr} -s -i stdin -g <(faidx {input.fasta} -i chromsizes) | \
        sort -k1,1 -k2,2n | \
        bedtools merge -s -d -1 -c 4,5,6,7,8,9 -o collapse,max,first,max,max,max -i stdin | \
        sort -k1,1 -k2,2n | \
        bedtools intersect -a stdin -b {input.motifs} -sorted -F 1 -wao | \
        cut -f18 --complement | \
        cat <(echo -e "chrom\tregion_start\tregion_end\tpeak_id\tpeak_score\tpeak_strand\tpeak_lfc\tpeak_logpval\tpeak_logqval\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | \
        pigz -f > {output}) &> {log}
        """

rule get_random_motifs:
    input:
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        motifs = build_annotations("motifs/" + config["genome"]["name"] + "_allmotifs.bed")
    output:
        "motifs/random_sequences-allFIMOresults.tsv.gz"
    params:
        window = config["motifs"]["enrichment-upstream"] + config["motifs"]["enrichment-downstream"] + 1,
        n = 6000
    log: "logs/get_random_motifs.log"
    shell: """
        (bedtools random -l {params.window} -n {params.n} -g <(faidx {input.fasta} -i chromsizes) | \
        sort -k1,1 -k2,2n | \
        bedtools merge -s -d -1 -c 4,5,6 -o collapse,first,first -i stdin | \
        sort -k1,1 -k2,2n | \
        bedtools intersect -a stdin -b {input.motifs} -sorted -F 1 -wao | \
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, 0, 0, 0, $7, $8, $9, $10, $11, $12, $13, $14}}' | \
        cat <(echo -e "chrom\tregion_start\tregion_end\tpeak_id\tpeak_score\tpeak_strand\tpeak_lfc\tpeak_logpval\tpeak_logqval\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | \
        pigz -f > {output}) &> {log}
        """

rule test_motif_enrichment:
    input:
        fimo_pos = "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-allFIMOresults.tsv.gz",
        fimo_neg = lambda wc: "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-unchanged-allFIMOresults.tsv.gz".format(**wc) if wc.negative=="unchanged" else "motifs/random_sequences-allFIMOresults.tsv.gz",
    output:
        tsv = "motifs/{condition}-v-{control}/{norm}/{category}/{negative}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-v-{negative}-motif_enrichment.tsv",
        plot = "motifs/{condition}-v-{control}/{norm}/{category}/{negative}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-v-{negative}-motif_enrichment.svg",
    params:
        fdr = config["motifs"]["enrichment-fdr"],
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated"
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/motif_enrichment.R"

#0. extend peak summit annotation to upstream and downstream distances
#1. if multiple annotations overlap on same strand, keep the one that is the most significant (avoid multiple-counting poorly called peaks erroneously split into multiple peaks)
rule get_meme_sequences:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}-summits.bed",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.fa"
    params:
        upstr = config["motifs"]["meme-chip"]["upstream"],
        dnstr = config["motifs"]["meme-chip"]["downstream"]
    log: "logs/get_meme_sequences/get_meme_sequences_{condition}-v-{control}-{norm}-{category}-{direction}.log"
    shell: """
        (bedtools slop -l {params.upstr} -r {params.dnstr} -s -i {input.peaks} -g <(faidx {input.fasta} -i chromsizes) | \
        bedtools cluster -s -d 0 -i stdin | \
        bedtools groupby -g 7 -c 5 -o max -full -i stdin | \
        sort -k4,4V | \
        bedtools getfasta -name+ -s -fi {input.fasta} -bed stdin > {output}) &> {log}
        """

rule meme_chip:
    input:
        seq = "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.fa",
        genome_fasta = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        dbs = build_annotations("motifs/" + config["genome"]["name"] + "_allmotifs.meme") if config["motifs"]["databases"] else []
    output:
        "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}-meme_chip/summary.tsv"
    params:
        db_command = "-db" if config["motifs"]["databases"] else [],
        meme_mode = config["motifs"]["meme-chip"]["meme-mode"],
        meme_nmotifs = config["motifs"]["meme-chip"]["meme-nmotifs"],
    log: "logs/meme_chip/meme_chip_{condition}-v-{control}-{norm}-{category}-{direction}.log"
    # threads: 4
    shell: """
        (meme-chip -oc motifs/{wildcards.condition}-v-{wildcards.control}/{wildcards.norm}/{wildcards.category}/{wildcards.condition}-v-{wildcards.control}_tss-seq-{wildcards.norm}-diffexp-results-{wildcards.category}-{wildcards.direction}-meme_chip {params.db_command} {input.dbs} -bfile <(fasta-get-markov {input.genome_fasta} -m 1) -order 1 -meme-mod {params.meme_mode} -meme-nmotifs {params.meme_nmotifs} -meme-p 1 -meme-norand -centrimo-local {input.seq}) &> {log}
        """

