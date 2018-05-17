#!/usr/bin/env python

#in this order: remove adapter, remove 3' molecular barcode, do NextSeq quality trimming
#reads shorter than 6nt after cleaning are discarded
rule clean_reads:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    output:
        fq = temp("fastq/cleaned/{sample}_tss-seq-trimmed.fastq.gz"),
        adapter = "logs/clean_reads/remove_adapter-{sample}.log",
        qual_trim = "logs/clean_reads/remove_3p_bc_and_trim-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"]
    threads: config["threads"]
    shell: """
        cutadapt --adapter={params.adapter} --error-rate=0.1 --minimum-length=18 --cores={threads} {input} 2> {output.adapter} | cutadapt --cut=-6 --nextseq-trim={params.trim_qual} --minimum-length=12 --output={output.fq} --cores={threads} - &> {output.qual_trim}
        """

rule remove_molecular_barcode:
    input:
        "fastq/cleaned/{sample}_tss-seq-trimmed.fastq.gz"
    output:
        fq = "fastq/cleaned/{sample}_tss-seq-clean.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    # threads: config["threads"]
    log: "logs/remove_molecular_barcode/remove_molecular_barcode-{sample}.log"
    shell: """
        (python scripts/extractMolecularBarcode.py {input} {output.fq} {output.barcodes} {output.ligation}) &> {log}
        """

