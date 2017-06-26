#!/usr/bin/env python
import os

configfile: "config.yaml"

CONTROL = config["samples"]["control"]
CONDITION = config["samples"]["condition"]

controlgroups = list(set([CONTROL[x]['group'] for x in CONTROL]))
conditiongroups = list(set([CONDITION[x]['group'] for x in CONDITION]))

#if using python >3.5:
#SAMPLES = {**CONTROL, **CONDITION}
SAMPLES = CONTROL.copy()
SAMPLES.update(CONDITION)
#SAMPLES = config["samples"]


localrules: all,
            make_stranded_genome,
            make_stranded_bedgraph,
            make_stranded_annotations,
            make_bigwig_for_deeptools,
            gzip_deeptools_table,
            gzip_deeptools_lfc_table,
            melt_matrix,
            melt_lfc_matrix,
            cat_matrices,
            cat_lfc_matrices,
            make_window_files,
            cat_windows,
            union_bedgraph,
            cat_strands

rule all:
    input:
        expand("qual_ctrl/fastqc/raw/{sample}", sample=SAMPLES),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip", sample=SAMPLES),
        expand("coverage/libsizenorm/{sample}-tss-libsizenorm-minus.bedgraph", sample=SAMPLES),
        expand("datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bygroup.png", annotation = config["annotations"], norm = ["spikenorm", "libsizenorm"], strand = ["SENSE", "ANTISENSE"]),
        #expand("coverage/{norm}/bw/lfc/{condition}-v-{control}-{norm}-{strand}.bw", norm=["spikenorm", "libsizenorm"], strand= ["SENSE", "ANTISENSE"], condition = CONDITION, control = CONTROL)
        #expand("datavis/{annotation}/{norm}/lfc/tss-lfc-{annotation}-{norm}-{strand}-heatmap-bygroup.png", annotation = config["annotations"], norm=["spikenorm", "libsizenorm"], strand = ["SENSE", "ANTISENSE"]),
        expand("correlations/{norm}-window{windowsize}-pca-scree.png", norm=["libsizenorm", "spikenorm"], windowsize=config["corr-binsizes"] ),
        #"diff_exp/de_bases/sig-bases-spikenorm.tsv",
        #"diff_exp/de_bases/sig-bases-libsizenorm.tsv"

rule fastqc_raw:
    input: 
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        "qual_ctrl/fastqc/raw/{sample}"
    threads: config["threads"]
    log: "logs/fastqc/raw/fastqc-raw-{sample}.log"
    shell: """
        #mkdir -p qual_ctrl/fastqc/raw/{wildcards.sample} 
        mkdir -p {output} 
        (fastqc -o {output} --noextract -t {threads} {input}) &> {log}
        """

#in this order: remove adapter, remove 3' molecular barcode, do NextSeq quality trimming
#reads shorter than 18 are thrown out, as the first 6 bases are the molecular barcode and 12-mer is around the theoretical minimum length to map uniquely to the combined Sc+Sp genome (~26 Mb)
rule remove_adapter:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        temp("fastq/cleaned/{sample}-noadapter.fastq.gz")
    params:
        adapter = config["cutadapt"]["adapter"],
    log: "logs/remove_adapter/remove_adapter-{sample}.log"
    shell: """
        (cutadapt -a {params.adapter} -m 24 -o {output} {input}) &> {log}
        """

rule remove_3p_barcode_and_qual_trim:
    input:
        "fastq/cleaned/{sample}-noadapter.fastq.gz"
    output:
        temp("fastq/cleaned/{sample}-trim.fastq")
    params:
        trim_qual = config["cutadapt"]["trim_qual"]
    log: "logs/remove_3p_bc_and_trim/cutadapt-{sample}.log"
    shell: """
        (cutadapt -u -6 --nextseq-trim={params.trim_qual} -m 18 -o {output} {input}) &> {log}
        """

rule remove_molec_barcode:
    input:
        "fastq/cleaned/{sample}-trim.fastq"
    output:
        fq = "fastq/cleaned/{sample}-clean.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    threads: config["threads"]
    log: "logs/remove_molec_barcode/removeMBC-{sample}.log"
    shell: """
        (python scripts/extractMolecularBarcode.py {input} fastq/cleaned/{wildcards.sample}-clean.fastq {output.barcodes} {output.ligation}) &> {log}
        pigz -f fastq/cleaned/{wildcards.sample}-clean.fastq
        """

rule fastqc_cleaned:
    input:
        "fastq/cleaned/{sample}-clean.fastq.gz"
    output:
        html = "qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.html",
        folder  = "qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip"
    threads : config["threads"]
    log: "logs/fastqc/cleaned/fastqc-cleaned-{sample}.log"
    shell: """
        mkdir -p qual_ctrl/fastqc/cleaned/{wildcards.sample}
        (fastqc -o qual_ctrl/fastqc/cleaned/{wildcards.sample} --noextract -t {threads} {input}) &> {log}
        """

rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"]
    output:
        expand("../genome/bowtie2_indexes/{basename}.{num}.bt2", basename=config["combinedgenome"]["name"], num = [1,2,3,4]),
        expand("../genome/bowtie2_indexes/{basename}.rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2])
    params:
        name = config["combinedgenome"]["name"]
    shell: """
        bowtie2-build {input.fasta} ../genome/bowtie2_indexes/{params.name} 
        """

rule align:
    input:
        expand("../genome/bowtie2_indexes/{basename}.{num}.bt2", basename=config["combinedgenome"]["name"], num = [1,2,3,4]),
        expand("../genome/bowtie2_indexes/{basename}.rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2]),
        gff = config["combinedgenome"]["gff"],
        fastq = "fastq/cleaned/{sample}-clean.fastq.gz"
    output:
        "alignment/{sample}/accepted_hits.bam"
    params:
        basename = config["combinedgenome"]["name"],
        read_mismatches = config["tophat2"]["read-mismatches"],
        read_gap_length = config["tophat2"]["read-gap-length"],
        read_edit_dist = config["tophat2"]["read-edit-dist"],
        min_anchor_length = config["tophat2"]["min-anchor-length"],
        splice_mismatches = config["tophat2"]["splice-mismatches"],
        min_intron_length = config["tophat2"]["min-intron-length"],
        max_intron_length = config["tophat2"]["max-intron-length"],
        max_insertion_length = config["tophat2"]["max-insertion-length"],
        max_deletion_length = config["tophat2"]["max-deletion-length"],
        max_multihits = config["tophat2"]["max-multihits"],
        segment_mismatches = config["tophat2"]["segment-mismatches"],
        segment_length = config["tophat2"]["segment-length"],
        min_coverage_intron = config["tophat2"]["min-coverage-intron"],
        max_coverage_intron = config["tophat2"]["max-coverage-intron"],
        min_segment_intron = config["tophat2"]["min-segment-intron"],
        max_segment_intron = config["tophat2"]["max-segment-intron"],
    conda:
        "envs/tophat2.yaml"
    threads : config["threads"]
    log: "logs/align/align-{sample}.log"
    shell: """
        (tophat2 --read-mismatches {params.read_mismatches} --read-gap-length {params.read_gap_length} --read-edit-dist {params.read_edit_dist} -o alignment/{wildcards.sample} --min-anchor-length {params.min_anchor_length} --splice-mismatches {params.splice_mismatches} --min-intron-length {params.min_intron_length} --max-intron-length {params.max_intron_length} --max-insertion-length {params.max_insertion_length} --max-deletion-length {params.max_deletion_length} --num-threads {threads} --max-multihits {params.max_multihits} --library-type fr-firststrand --segment-mismatches {params.segment_mismatches} --no-coverage-search --segment-length {params.segment_length} --min-coverage-intron {params.min_coverage_intron} --max-coverage-intron {params.max_coverage_intron} --min-segment-intron {params.min_segment_intron} --max-segment-intron {params.max_segment_intron} --b2-sensitive -G {input.gff} ../genome/bowtie2_indexes/{params.basename} {input.fastq}) &> {log}
        """

rule select_unique_mappers:
    input:
        "alignment/{sample}/accepted_hits.bam"
    output:
        temp("alignment/{sample}-unique.bam")
    threads: config["threads"]
    log: "logs/select_unique_mappers/select_unique_mappers-{sample}.log"
    shell: """
        (samtools view -b -h -q 50 -@ {threads} {input} | samtools sort -@ {threads} - > {output}) &> {log}
        """   

rule remove_PCR_duplicates:
    input:
        "alignment/{sample}-unique.bam"
    output:
        "alignment/{sample}-noPCRdup.bam"
    log: "logs/remove_PCR_duplicates/removePCRduplicates-{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output}) &> {log}
        """

rule get_coverage:
    input:
        "alignment/{sample}-noPCRdup.bam"
    output:
        SIplmin = "coverage/counts/spikein/{sample}-tss-SI-counts-plmin.bedgraph",
        SIpl = "coverage/counts/spikein/{sample}-tss-SI-counts-plus.bedgraph",
        SImin = "coverage/counts/spikein/{sample}-tss-SI-counts-minus.bedgraph",
        plmin = "coverage/counts/{sample}-tss-counts-plmin.bedgraph",
        plus = "coverage/counts/{sample}-tss-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-tss-counts-minus.bedgraph"
    params:
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
        si_prefix = config["combinedgenome"]["spikein_prefix"]
    log: "logs/get_coverage/get_coverage-{sample}.log"
    shell: """
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SIplmin}) &> {log};
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SIpl}) &>> {log};
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SImin}) &>> {log};
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.plmin}) &>> {log};
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.plus}) &>> {log};
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.minus}) &>> {log};
        """

rule normalize:
    input:
        plus = "coverage/counts/{sample}-tss-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-tss-counts-minus.bedgraph",
        plmin = "coverage/counts/{sample}-tss-counts-plmin.bedgraph",
        SIplmin = "coverage/counts/spikein/{sample}-tss-SI-counts-plmin.bedgraph"
    output:
        spikePlus = "coverage/spikenorm/{sample}-tss-spikenorm-plus.bedgraph",
        spikeMinus = "coverage/spikenorm/{sample}-tss-spikenorm-minus.bedgraph",
        libnormPlus = "coverage/libsizenorm/{sample}-tss-libsizenorm-plus.bedgraph",
        libnormMinus = "coverage/libsizenorm/{sample}-tss-libsizenorm-minus.bedgraph"
    log: "logs/normalize/normalize-{sample}.log"
    shell: """
        (scripts/libsizenorm.awk {input.SIplmin} {input.plus} > {output.spikePlus}) &> {log} 
        (scripts/libsizenorm.awk {input.SIplmin} {input.minus} > {output.spikeMinus}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.plus} > {output.libnormPlus}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.minus} > {output.libnormMinus}) &>> {log}
        """

#make 'stranded' genome for datavis purposes
rule make_stranded_genome:
    input:
        config["genome"]["chrsizes"]
    output:
        os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"
    log: "logs/make_stranded_genome.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input} > {output}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = "coverage/{norm}/{sample}-tss-{norm}-plus.bedgraph",        
        minus = "coverage/{norm}/{sample}-tss-{norm}-minus.bedgraph"        
    output:
        sense = "coverage/{norm}/{sample}-tss-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}-tss-{norm}-ANTISENSE.bedgraph"
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2, $3, $4}}' {input.plus} > coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-plus.tmp) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-minus", $2, $3, $4}}' {input.minus} > coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-minus.tmp) &>> {log}
        (cat coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-plus.tmp coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-minus.tmp | LC_COLLATE=C sort -k1,1 -k2,2n > {output.sense}) &>> {log} 
        rm coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-*.tmp 
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2, $3, $4}}' {input.minus} > coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-plus.tmp) &>> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-minus", $2, $3, $4}}' {input.plus} > coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-minus.tmp) &>> {log}
        (cat coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-plus.tmp coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-minus.tmp | LC_COLLATE=C sort -k1,1 -k2,2n > {output.antisense}) &>> {log} 
        rm coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-*.tmp 
        """

rule make_stranded_annotations:
    input:
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "../genome/annotations/stranded/{annotation}-STRANDED.bed"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}$6=="+"{{print $1"-plus", $2, $3, $4, $5, $6}} $6=="-"{{print $1"-minus", $2, $3, $4, $5, $6}}' {input} > {output}) &> {log}
        """

rule make_bigwig_for_deeptools:
    input:
        bedgraph = "coverage/{norm}/{sample}-tss-{norm}-{strand}.bedgraph",
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"
    output:
        "coverage/{norm}/bw/{sample}-tss-{norm}-{strand}.bw",
    log : "logs/make_bigwig_for_deeptools/make_bigwig_for_deeptools-{sample}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

rule bigwig_compare:
    input:
        condition = "coverage/{norm}/bw/{condition}-tss-{norm}-{strand}.bw",
        control = "coverage/{norm}/bw/{control}-tss-{norm}-{strand}.bw"
    output:
        "coverage/{norm}/bw/lfc/{condition}-v-{control}-{norm}-{strand}.bw"
    log: "logs/bigwig_compare/bigwig_compare-{condition}-v-{control}-{norm}-{strand}.log"
    threads: config["threads"]
    shell: """
        (bigwigCompare -b1 {input.condition} -b2 {input.control} --pseudocount 0.1 --ratio log2 --binSize 1 -p {threads} -o {output}) &> {log}
        """

rule deeptools_matrix:
    input:
        annotation = "../genome/annotations/stranded/{annotation}-STRANDED.bed",
        bw = "coverage/{norm}/bw/{sample}-tss-{norm}-{strand}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}-{norm}-{strand}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}
        """

rule deeptools_lfc_matrix:
    input:
        annotation = "../genome/annotations/stranded/{annotation}-STRANDED.bed",
        bw = "coverage/{norm}/bw/lfc/{condition}-v-{control}-{norm}-{strand}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/lfc/{annotation}-{condition}-v-{control}-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/lfc/{annotation}-{condition}-v-{control}-{norm}-{strand}.tsv")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/compute_lfc_Matrix-{annotation}-{condition}-v-{control}-{norm}-{strand}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}
        """

rule gzip_deeptools_table:
    input:
        tsv = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv",
        mat = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.mat.gz"
    output:
        "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv.gz"
    shell: """
        pigz -f {input.tsv}
        rm {input.mat}
        """

rule gzip_deeptools_lfc_table:
    input:
        tsv = "datavis/{annotation}/{norm}/lfc/{annotation}-{condition}-v-{control}-{norm}-{strand}.tsv",
        mat = "datavis/{annotation}/{norm}/lfc/{annotation}-{condition}-v-{control}-{norm}-{strand}.mat.gz"
    output:
        "datavis/{annotation}/{norm}/lfc/{annotation}-{condition}-v-{control}-{norm}-{strand}.tsv.gz"
    shell: """
        pigz -f {input.tsv}
        rm {input.mat}
        """

rule melt_matrix:
    input:
        matrix = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv.gz"
    output:
        temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}-melted.tsv.gz")
    params:
        name = lambda wildcards : wildcards.sample,
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"]
    script:
        "scripts/melt_matrix2.R"

rule melt_lfc_matrix:
    input:
        matrix = "datavis/{annotation}/{norm}/lfc/{annotation}-{condition}-v-{control}-{norm}-{strand}.tsv.gz"
    output:
        "datavis/{annotation}/{norm}/lfc/{annotation}-{condition}-v-{control}-{norm}-{strand}-melted.tsv.gz"
    params:
        condition = lambda wildcards : wildcards.condition,
        control = lambda wildcards : wildcards.control,
        conditiongroup = lambda wildcards : SAMPLES[wildcards.condition]["group"],
        controlgroup = lambda wildcards : SAMPLES[wildcards.control]["group"],
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"]
    script:
        "scripts/melt_lfc_matrix.R"

rule cat_matrices:
    input:
        expand("datavis/{{annotation}}/{{norm}}/{{annotation}}-{sample}-{{norm}}-{{strand}}-melted.tsv.gz", sample=SAMPLES)

    output:
        "datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{strand}.tsv.gz"
    shell: """
        cat {input} > {output}
        """

rule cat_lfc_matrices:
    input:
        expand("datavis/{{annotation}}/{{norm}}/lfc/{{annotation}}-{condition}-v-{control}-{{norm}}-{{strand}}-melted.tsv.gz", condition = CONDITION, control = CONTROL)
    output:
        "datavis/{annotation}/{norm}/lfc/allsampleslfc-{annotation}-{norm}-{strand}.tsv.gz"
    shell: """
        cat {input} > {output}
        """

rule r_datavis:
    input:
        matrix = "datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{strand}.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bysample.png",
        heatmap_group = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bygroup.png",
    params:
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        #pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
        #metagene_palette = lambda wildcards : config["annotations"][wildcards.annotation]["metagene_palette"],
        #avg_heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["avg_heatmap_cmap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plotHeatmaps.R"


rule r_lfc_datavis:
    input:
        matrix = "datavis/{annotation}/{norm}/lfc/allsampleslfc-{annotation}-{norm}-{strand}.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/{norm}/lfc/tss-lfc-{annotation}-{norm}-{strand}-heatmap-bysample.png",
        heatmap_group = "datavis/{annotation}/{norm}/lfc/tss-lfc-{annotation}-{norm}-{strand}-heatmap-bygroup.png",
    params:
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        #pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["lfc_heatmap_colormap"],
        #metagene_palette = lambda wildcards : config["annotations"][wildcards.annotation]["metagene_palette"],
        #avg_heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["avg_heatmap_cmap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plot_lfc_Heatmaps.R"

rule make_window_files:
    input:
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"
    output:
        temp(os.path.dirname(config["genome"]["chrsizes"]) + "/windows-{windowsize}.bed")
    log: "logs/make_window_files/make_window_files-{windowsize}.log"
    shell: """
        (bedtools makewindows -g {input.chrsizes} -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule map_to_windows:
    input:
        bed = os.path.dirname(config["genome"]["chrsizes"]) + "/windows-{windowsize}.bed",
        bedgraph = "coverage/{norm}/{sample}-tss-{norm}-SENSE.bedgraph"
    output:
        temp("coverage/{norm}/.{sample}-{norm}-{windowsize}.tsv")
    log: "logs/map_to_windows/map_to_windows-{sample}-{norm}-{windowsize}.log"
    shell: """
        (bedtools map -c 4 -o sum -a {input.bed} -b {input.bedgraph} | cut -f4 > {output}) &> {log}
        """
 
rule cat_windows:
    input:
        values = expand("coverage/{{norm}}/.{sample}-{{norm}}-{{windowsize}}.tsv", sample=SAMPLES),
        coord = os.path.dirname(config["genome"]["chrsizes"]) + "/windows-{windowsize}.bed",
    output:
        "correlations/{norm}-window{windowsize}.tsv"
    params:
        labels = list(SAMPLES.keys())
    shell: """
        echo -e "chr\tstart\tend\t{params.labels}\n$(paste {input.coord} {input.values})" > {output}
        """

rule plot_correlations:
    input:
        "correlations/{norm}-window{windowsize}.tsv"
    output:
        scatter = "correlations/{norm}-window{windowsize}-pairwise-scatter.png",
        dists_cluster = "correlations/{norm}-window{windowsize}-sample-dists-clustered.png",
        dists_nocluster = "correlations/{norm}-window{windowsize}-sample-dists-unclustered.png",
        pca = "correlations/{norm}-window{windowsize}-pca.png",
        scree = "correlations/{norm}-window{windowsize}-pca-scree.png"
    script:
        "scripts/plotcorrelations.R"

name_string = " ".join(SAMPLES)

rule union_bedgraph:
    input:
        exp = expand("coverage/counts/{sample}-tss-counts-{{strand}}.bedgraph", sample=SAMPLES),
        si = expand("coverage/counts/spikein/{sample}-tss-SI-counts-{{strand}}.bedgraph", sample=SAMPLES)
    output:
        exp = temp("coverage/counts/union-bedgraph-{strand}-nozero.txt"),    
        si = temp("coverage/counts/spikein/union-bedgraph-si-{strand}-nozero.txt")
    shell: """
        bedtools unionbedg -i {input.exp} -header -names {name_string} |
        awk -v awkstrand={wildcards.strand} 'BEGIN{{FS=OFS="\t"}}{{print awkstrand, $0 }}' |
        awk 'BEGIN{{FS=OFS="\t"}}{{t=0; for(i=5; i<=NF; i++) t+=$i}} t>1{{print $0}}' > {output.exp}
        bedtools unionbedg -i {input.si} -header -names {name_string} |
        awk -v awkstrand={wildcards.strand} 'BEGIN{{FS=OFS="\t"}}{{print awkstrand, $0 }}' |
        awk 'BEGIN{{FS=OFS="\t"}}{{t=0; for(i=5; i<=NF; i++) t+=$i}} t>1{{print $0}}' > {output.si}
        """

rule cat_strands:
    input:
        exp = expand("coverage/counts/union-bedgraph-{strand}-nozero.txt", strand=["plus", "minus"]),
        si = expand("coverage/counts/spikein/union-bedgraph-si-{strand}-nozero.txt", strand=["plus", "minus"])
    output:
        exp = "coverage/counts/union-bedgraph-bothstr-nozero.txt",
        si = "coverage/counts/spikein/union-bedgraph-si-bothstr-nozero.txt"
    shell: """
        cat {input.exp} > coverage/counts/.catstrandtemp.txt
        cut -f1-4 coverage/counts/.catstrandtemp.txt | awk 'BEGIN{{FS="\t"; OFS=":"}}{{print $1, $2, $3, $4}}'  > .positions.txt
        cut -f5- coverage/counts/.catstrandtemp.txt > .values.txt
        paste .positions.txt .values.txt > {output.exp}
        rm coverage/counts/.catstrandtemp.txt .positions.txt .values.txt
        cat {input.si} > coverage/counts/.sicatstrandtemp.txt
        cut -f1-4 coverage/counts/.sicatstrandtemp.txt | awk 'BEGIN{{FS="\t"; OFS=":"}}{{print $1, $2, $3, $4}}' > .sipositions.txt
        cut -f5- coverage/counts/.sicatstrandtemp.txt > .sivalues.txt
        paste .sipositions.txt .sivalues.txt > {output.si}
        rm coverage/counts/.sicatstrandtemp.txt .sipositions.txt .sivalues.txt
        """

rule get_de_bases:
   input:
        exp = "coverage/counts/union-bedgraph-bothstr-nozero.txt",
        si = "coverage/counts/spikein/union-bedgraph-si-bothstr-nozero.txt"
   params:
        alpha = config["deseq"]["fdr"],
        samples = list(SAMPLES.keys()),
        samplegroups = [SAMPLES[x]["group"] for x in SAMPLES]
   output:
        spikenorm = protected("diff_exp/de_bases/sig-bases-spikenorm.tsv"),
        libsizenorm = protected("diff_exp/de_bases/sig-bases-libsizenorm.tsv")
   script:
        "scripts/base-diff-expr.R"



