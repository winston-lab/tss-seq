#!/usr/bin/env python
import os.path

configfile: "config.yaml"

CONTROL = config["control"]
CONDITION = config["condition"]

#if using python >3.5:
#SAMPLES = {**CONTROL, **CONDITION}
SAMPLES = CONTROL.copy()
SAMPLES.update(config["condition"])

FNAMES = list(map(lambda path : os.path.basename(path).split('.')[0] , list(SAMPLES.values())))



localrules: all,

rule all:
    input:
        #expand("qual_ctrl/fastqc/raw/{sample}/{fname}_fastqc.zip", sample=SAMPLES, fname = SAMPLES[sample]),
        expand("qual_ctrl/fastqc/raw/{sample}", sample=SAMPLES),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip", sample=SAMPLES),
        expand("coverage/libsizenorm/{sample}.libsizenorm.minus.bedgraph", sample=SAMPLES)

rule fastqc_raw:
    input: 
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        #html = "qual_ctrl/fastqc/raw/{sample}/{fname}_fastqc.html",
        #folder = "qual_ctrl/fastqc/raw/{sample}/{fname}_fastqc.zip",
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
        lambda wildcards: SAMPLES[wildcards.sample]
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
        SCplmin = "coverage/counts/scer/{sample}.SC.counts.plmin.bedgraph",
        SCpl = "coverage/counts/scer/{sample}.SC.counts.plus.bedgraph",
        SCmin = "coverage/counts/scer/{sample}.SC.counts.minus.bedgraph",
        SPplmin = "coverage/counts/{sample}.counts.plmin.bedgraph",
        SPpl = "coverage/counts/{sample}.counts.plus.bedgraph",
        SPmin = "coverage/counts/{sample}.counts.minus.bedgraph"
    log: "logs/get_coverage/get_coverage-{sample}.log"
    shell: """
        (genomeCoverageBed -bga -5 -ibam {input} | grep Scer_ | sed 's/Scer_//g' | sort -k1,1 -k2,2n > {output.SCplmin}) &> {log};
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep Scer_ | sed 's/Scer_//g' | sort -k1,1 -k2,2n > {output.SCpl}) &>> {log};
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep Scer_ | sed 's/Scer_//g' | sort -k1,1 -k2,2n > {output.SCmin}) &>> {log};
        (genomeCoverageBed -bga -5 -ibam {input} | grep Spom_ | sed 's/Spom_//g' | sort -k1,1 -k2,2n > {output.SPplmin}) &>> {log};
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep Spom_ | sed 's/Spom_//g' | sort -k1,1 -k2,2n > {output.SPpl}) &>> {log};
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep Spom_ | sed 's/Spom_//g' | sort -k1,1 -k2,2n > {output.SPmin}) &>> {log};
        """

rule normalize:
    input:
        SPpl = "coverage/counts/{sample}.counts.plus.bedgraph",
        SPmin = "coverage/counts/{sample}.counts.minus.bedgraph",
        SPplmin = "coverage/counts/{sample}.counts.plmin.bedgraph",
        SCplmin = "coverage/counts/scer/{sample}.SC.counts.plmin.bedgraph"
    output:
        spikePlus = "coverage/spikenorm/{sample}.spikenorm.plus.bedgraph",
        spikeMinus = "coverage/spikenorm/{sample}.spikenorm.minus.bedgraph",
        libnormPlus = "coverage/libsizenorm/{sample}.libsizenorm.plus.bedgraph",
        libnormMinus = "coverage/libsizenorm/{sample}.libsizenorm.minus.bedgraph"
    log: "logs/normalize/normalize-{sample}.log"
    shell: """
        (scripts/libsizenorm.awk {input.SCplmin} {input.SPpl} > {output.spikePlus}) &> {log} 
        (scripts/libsizenorm.awk {input.SCplmin} {input.SPmin} > {output.spikeMinus}) &>> {log}
        (scripts/libsizenorm.awk {input.SPplmin} {input.SPpl} > {output.libnormPlus}) &>> {log}
        (scripts/libsizenorm.awk {input.SPplmin} {input.SPmin} > {output.libnormMinus}) &>> {log}
        """

