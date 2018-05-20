#!/usr/bin/env python

#run fimo in parallel for each motif
# TODO: for maximum efficiency should chunk the motif list into
# something bigger than a single motif to avoid eating the
# overhead of waiting for cluster job submission
rule fimo:
    input:
        fasta = config["genome"]["fasta"],
        motif_db = config["motifs"]["databases"]
    params:
        alpha = config["motifs"]["fimo-pval"]
    output:
        tsv = temp("motifs/.{motif}.tsv"),
        bed = temp("motifs/.{motif}.bed") #first 6 columns are BED6, plus extra info in later columns
    shell: """
        fimo --motif {wildcards.motif} --bgfile <(fasta-get-markov {input.fasta}) --parse-genomic-coord --thresh {params.alpha} --text <(meme2meme -bg <(fasta-get-markov {input.fasta}) {input.motif_db}) {input.fasta} | sed -e 's/\//_/g; s/&/_/g; s/{{/[/g; s/}}/]/g' | tee {output.tsv} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $3, $4, $5+1, $1, -log($8)/log(10), $6, $2, $10}}' | sort -k1,1 -k2,2n > {output.bed}
        """

rule cat_fimo_motifs:
    input:
        tsv = expand("motifs/.{motif}.tsv", motif=MOTIFS),
        bed = expand("motifs/.{motif}.bed", motif=MOTIFS)
    output:
        tsv = "motifs/allmotifs.tsv.gz",
        bed = "motifs/allmotifs.bed",
    threads: config["threads"]
    shell: """
        cat {input.tsv} | pigz -f > {output.tsv}
        cat {input.bed} | sort -k1,1 -k2,2n --parallel={threads} > {output.bed}
        """

#bedtools intersect peaks with fimo motifs
#0. with the summit of the peak as reference, extend annotation to upstream and 'downstream' distances
#1. if multiple annotations overlap on same strand, keep the one that is the most significant (avoid multiple-counting, so that poorly called peaks which are erroneously split into multiple peaks do not bias the motif enrichment analyses.)
#2. intersect with motif file
rule get_upstream_motifs:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.narrowpeak",
        chrsizes = config["genome"]["chrsizes"],
        motifs = "motifs/allmotifs.bed"
    output:
        "motifs/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-motifs.tsv.gz"
    params:
        upstr = config["motifs"]["enrichment-upstream"],
        dnstr = config["motifs"]["enrichment-downstream"]
    log: "logs/get_upstream_motifs/get_upstream_motifs-{condition}-v-{control}-{norm}-{category}-{direction}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' {input.peaks} | bedtools slop -l {params.upstr} -r {params.dnstr} -s -i stdin -g {input.chrsizes} | bedtools cluster -s -d 0 -i stdin | bedtools groupby -g 11 -c 9 -o max -full -i stdin | cut -f1-10 | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b {input.motifs} -sorted -wao | cut -f19 --complement | cat <(echo -e "chrom\tregion_start\tregion_end\tpeak_id\tpeak_score\tpeak_strand\tpeak_lfc\tpeak_logpval\tpeak_logqval\tpeak_summit\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | pigz -f > {output}) &> {log}
        """

rule get_random_motifs:
    input:
        chrsizes = config["genome"]["chrsizes"],
        motifs = "motifs/allmotifs.bed"
    output:
        "motifs/random-motifs.tsv.gz"
    params:
        n = 6000,
        size = config["motifs"]["enrichment-upstream"] + config["motifs"]["enrichment-downstream"]
    log: "logs/get_random_motifs.log"
    shell: """
        (bedtools random -l {params.size} -n {params.n} -g {input.chrsizes} | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b {input.motifs} -sorted -wao | cut -f15 --complement | cat <(echo -e "chrom\tregion_start\tregion_end\tregion_id\tregion_length\tregion_strand\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | pigz -f > {output}) &> {log}
        """

rule test_motif_enrichment:
    input:
        fimo_pos = "motifs/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-motifs.tsv.gz",
        fimo_neg = lambda wc: "motifs/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-{category}-unchanged-motifs.tsv.gz".format(**wc) if wc.negative=="unchanged" else "motifs/random-motifs.tsv.gz"
    params:
        pval_cutoff = config["motifs"]["fimo-pval"],
        alpha= config["motifs"]["enrichment-fdr"],
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated"
    output:
        tsv = "motifs/{condition}-v-{control}/{norm}/{negative}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-v-{negative}-motif_enrichment.tsv",
        plot = "motifs/{condition}-v-{control}/{norm}/{negative}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-v-{negative}-motif_enrichment.svg",
    script: "../scripts/motif_enrichment.R"

#NOTE: below are rules for visualizing motif occurrences, but need to have a way to do this efficiently/in an interpretable way for thousands of motifs
# rule get_motif_coverage:
#     input:
#         bed = "motifs/.{motif}.bed", #this is sorted when created
#         chrsizes = config["genome"]["chrsizes"]
#     output:
#         bg = "motifs/coverage/{motif}.bedgraph",
#         bw = "motifs/coverage/{motif}.bw",
#     shell: """
#         cut -f1-6 {input.bed} | bedtools genomecov -bga -i stdin -g {input.chrsizes} | sort -k1,1 -k2,2n > {output.bg}
#         bedGraphToBigWig {output.bg} {input.chrsizes} {output.bw}
#         """

# rule motif_matrix:
#     input:
#         annotation = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed",
#         bw = "motifs/coverage/{motif}.bw"
#     output:
#         dtfile = temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.mat"),
#         matrix = temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv"),
#         matrix_gz = "motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv.gz",
#     params:
#         refpoint = "TSS",
#         upstream = config["motifs"]["upstream"] + config["motifs"]["binsize"],
#         dnstream = config["motifs"]["freq-downstream"] + config["motifs"]["binsize"],
#         binsize = config["motifs"]["binsize"],
#         sort = "keep",
#         binstat = "sum"
#     threads : config["threads"]
#     log: "logs/deeptools/computeMatrix-{motif}-{condition}-{control}-{norm}-{direction}-{category}.log"
#     shell: """
#         (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --averageTypeBins {params.binstat} -p {threads}) &> {log}
#         pigz -fk {output.matrix}
#         """

# rule melt_motif_matrix:
#     input:
#         matrix = "motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv.gz",
#     output:
#         temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}-melted.tsv.gz"),
#     params:
#         refpoint = "TSS",
#         binsize = config["motifs"]["binsize"],
#         upstream = config["motifs"]["upstream"],
#     script:
#         "scripts/melt_motif_matrix.R"

# rule cat_motif_matrices:
#     input:
#         expand("motifs/datavis/{motif}_{{condition}}-v-{{control}}_{{norm}}-{direction}-peaks-{category}-melted.tsv.gz", category=CATEGORIES, motif=MOTIFS, direction=["up","down","unchanged"]),
#     output:
#         "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.tsv.gz"
#     log: "logs/cat_matrices/cat_matrices-{condition}-{control}-{norm}.log"
#     shell: """
#         (cat {input} > {output}) &> {log}
#         """

rule plot_motif_freq:
    input:
        "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.tsv.gz"
    output:
        "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.svg"
    script: "scripts/motif_metagenes.R"

#0. filter out double counted peaks (sometimes a peak can be 'genic' for two genes, causing it to be listed twice)
#1. for sequences upstream of peaks: with the START of the peak as reference, extend annotation to upstream and 'downstream' distances; for peak sequences, just take the peak sequence
#2. for sequences upstream of peaks: if multiple annotations overlap on same strand, keep the one that is the most significant (avoid multiple-counting poorly called peaks erroneously split into multiple peaks); for peak sequences, no such limitation since they should be non-overlapping
rule get_meme_de_peak_sequences:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed",
        chrsizes = config["genome"]["chrsizes"],
        fasta = config["genome"]["fasta"]
    output:
        "motifs/meme/{condition}-v-{control}/{norm}/{region}/{condition}-v-{control}-results-{norm}-{direction}-{category}-{region}-meme.fa"
    params:
        upstr = config["motifs"]["meme-chip"]["upstream"],
        dnstr = config["motifs"]["meme-chip"]["downstream"]
    log: "logs/get_meme_de_peak_sequences/get_meme_de_peak_sequences-{condition}-v-{control}-{norm}-{direction}-{category}-{region}.log"
    run:
        if wildcards.region=="upstream":
            shell("""(uniq {input.peaks} | bedtools flank -l {params.upstr} -r 0 -s -i stdin -g {input.chrsizes} | bedtools slop -l 0 -r {params.dnstr} -s -i stdin -g {input.chrsizes} | bedtools cluster -s -d 0 -i stdin | sed 's/:/\t/g' | bedtools groupby -g 8 -c 6 -o max -full -i stdin | sort -k6,6nr | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5":"$6, $7}}' | bedtools getfasta -name+ -s -fi {input.fasta} -bed stdin > {output}) &> {log}""")
        elif wildcards.region=="peak":
            shell("""(uniq {input.peaks} | bedtools getfasta -name+ -s -fi {input.fasta} -bed stdin > {output}) &> {log}""")

rule meme_chip:
    input:
        seq = "motifs/meme/{condition}-v-{control}/{norm}/{region}/{condition}-v-{control}-results-{norm}-{direction}-{category}-{region}-meme.fa",
        genome_fasta = config["genome"]["fasta"],
        dbs = config["motifs"]["databases"]
    output:
        "motifs/meme/{condition}-v-{control}/{norm}/{region}/{condition}-v-{control}-results-{norm}-{direction}-{category}-{region}-meme_chip/meme-chip.html"
    params:
        dbs = ["-db " + x for x in config["motifs"]["databases"]],
        nmeme = lambda wc: int(1e5//(config["motifs"]["meme-chip"]["upstream"] + config["motifs"]["meme-chip"]["downstream"])) if wc.region=="upstream" else int(1e5//config["motifs"]["meme-chip"]["peak-ccut"]),
        ccut = lambda wc: int(config["motifs"]["meme-chip"]["upstream"] + config["motifs"]["meme-chip"]["downstream"]) if wc.region=="upstream" else config["motifs"]["meme-chip"]["peak-ccut"],
        meme_mode = config["motifs"]["meme-chip"]["meme-mode"],
        meme_nmotifs = config["motifs"]["meme-chip"]["meme-nmotifs"],
    conda: "envs/meme_chip.yaml"
    shell: """
        meme-chip -oc motifs/meme/{wildcards.condition}-v-{wildcards.control}/{wildcards.norm}/{wildcards.region}/{wildcards.condition}-v-{wildcards.control}-results-{wildcards.norm}-{wildcards.direction}-{wildcards.category}-{wildcards.region}-meme_chip -bfile <(fasta-get-markov {input.genome_fasta} -m 1) -nmeme {params.nmeme} -norand -ccut {params.ccut} -meme-mod {params.meme_mode} -meme-nmotifs {params.meme_nmotifs} -centrimo-local {params.dbs} {input.seq}
        """

