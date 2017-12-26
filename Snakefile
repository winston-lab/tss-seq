#!/usr/bin/env python
import os
import subprocess
from math import log2, log10

configfile: "config.yaml"

SAMPLES = config["samples"]
sisamples = {k:v for (k,v) in SAMPLES.items() if v["spikein"]=="y"}
PASSING = {k:v for (k,v) in SAMPLES.items() if v["pass-qc"] == "pass"}
sipassing = {k:v for (k,v) in PASSING.items() if v["spikein"] == "y"}

#groups which have at least two passing samples, so that they are valid for peakcalling and diff exp
validgroups = set([z for z in [PASSING[x]['group'] for x in PASSING] if [PASSING[x]['group'] for x in PASSING].count(z)>=2])
validgroups_si = set([z for z in [PASSING[x]['group'] for x in PASSING if PASSING[x]['spikein']=="y"] if [PASSING[x]['group'] for x in PASSING].count(z)>=2])
controlgroups = [g for g in config["comparisons"]["libsizenorm"]["controls"] if g in validgroups]
conditiongroups = [g for g in config["comparisons"]["libsizenorm"]["conditions"] if g in validgroups]
controlgroups_si = [g for g in config["comparisons"]["spikenorm"]["controls"] if g in validgroups_si]
conditiongroups_si = [g for g in config["comparisons"]["spikenorm"]["conditions"] if g in validgroups_si]

CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

localrules:
    all,
    bowtie2_build,
    get_si_pct, cat_si_pct, plot_si_pct,
    make_stranded_genome, make_stranded_annotations,
    tss_peaks_to_narrowpeak,
    build_genic_annotation, build_convergent_annotation, build_divergent_annotation, build_intergenic_annotation,
    classify_peaks_genic, classify_peaks_intragenic, classify_peaks_antisense,
    classify_peaks_convergent, classify_peaks_divergent, classify_peaks_intergenic,
    combine_tss_peaks, map_counts_to_peaks, get_peak_counts,
    separate_de_peaks, de_peaks_to_bed,
    get_de_genic, get_de_intragenic, get_de_antisense,
    get_de_convergent, get_de_divergent, get_de_intergenic,
    # get_de_intragenic_frequency
    # plot_de_intragenic_frequency
    # get_intra_orfs
    separate_sig_de, get_de_category_bed,
    get_peak_sequences_all,
    # meme_chip
    # class_v_genic

rule all:
    input:
        #FastQC
        expand("qual_ctrl/fastqc/raw/{sample}", sample=SAMPLES),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip", sample=SAMPLES),
        #alignment
        expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}-tss-{norm}-{strand}.bw", norm=["spikenorm","libsizenorm", "counts", "sicounts"], sample=SAMPLES, strand=["SENSE","ANTISENSE","plus","minus"]),
        #datavis
        expand(expand("datavis/{{annotation}}/libsizenorm/tss-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}-{{strand}}-heatmap-bygroup.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], status=["all","passing"], strand=["SENSE","ANTISENSE"]),
        expand(expand("datavis/{{annotation}}/spikenorm/tss-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}-{{strand}}-heatmap-bygroup.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), annotation=config["annotations"], status=["all","passing"], strand=["SENSE","ANTISENSE"]),
        #quality control
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all", "passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-tss-{{status}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status = ["all", "passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-tss-{{status}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status = ["all", "passing"]),
        #find intragenic ORFs
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intragenic-orfs/{condition}-v-{control}-libsizenorm-{{direction}}-intragenic-orfs.tsv", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"]),
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intragenic-orfs/{condition}-v-{control}-spikenorm-{{direction}}-intragenic-orfs.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"]),
        #peakcalling on all samples
        expand("peakcalling/{sample}-exp-allpeaks.narrowPeak", sample=SAMPLES),
        expand("peakcalling/{sample}-si-allpeaks.narrowPeak", sample=sisamples),
        #IDR for all groups which have at least two passing samples
        expand("peakcalling/{group}-exp-idrpeaks.narrowPeak", group=validgroups),
        expand("peakcalling/{group}-si-idrpeaks.narrowPeak", group=validgroups_si),
        #classify peaks into categories
        expand("peakcalling/{category}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
        expand("peakcalling/{condition}-v-{control}-peakdistances.svg", zip, condition=conditiongroups + ["all"], control=controlgroups + ["all"]),
        #combine called peaks for conditions vs control
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-peaks.bed", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-si-peaks.bed", zip, condition=conditiongroups_si, control=controlgroups_si),
        #differential expression of peaks
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-peak-counts.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-si-peak-counts.tsv", zip, condition=conditiongroups_si, control=controlgroups_si),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-libsizenorm-all.bed", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-spikenorm-all.bed", zip, condition=conditiongroups_si, control=controlgroups_si),
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-libsizenorm-{{dir}}.{{fmt}}", zip, condition=conditiongroups, control=controlgroups), dir=["up","unchanged","down"], fmt=["tsv","bed"]),
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-spikenorm-{{dir}}.{{fmt}}", zip, condition=conditiongroups_si, control=controlgroups_si), dir=["up","unchanged", "down"], fmt=["tsv","bed"]),
        #categorize DE peaks
        expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-results-libsizenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups, control=controlgroups), direction = ["all","up","unchanged","down"], category=CATEGORIES),
        expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-results-spikenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["all","up","unchanged","down"], category=CATEGORIES),
        #DE summary
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-libsizenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups, control=controlgroups), plot = ["summary", "maplot", "volcano"]),
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-spikenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), plot = ["summary", "maplot", "volcano"]),
        #intragenic frequency per ORF
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-libsizenorm-{{direction}}-freqperORF.svg", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"]),
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-spikenorm-{{direction}}-freqperORF.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"]),
        expand("diff_exp/{condition}-v-{control}/genic_v_class/{condition}-v-{control}-libsizenorm-genic-v-class.svg", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/genic_v_class/{condition}-v-{control}-spikenorm-genic-v-class.svg", zip, condition=conditiongroups_si, control=controlgroups_si),
        #FIMO
        "motifs/allmotifs.bed",
        # expand(expand("motifs/{condition}-v-{control}/{condition}-v-{control}_libsizenorm-{{direction}}-{{category}}-motifs.tsv.gz", zip, condition=conditiongroups, control=controlgroups), direction=["up","unchanged","down"], category=CATEGORIES),
        # expand(expand("motifs/{condition}-v-{control}/{condition}-v-{control}_spikenorm-{{direction}}-{{category}}-motifs.tsv.gz", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","unchanged","down"], category=CATEGORIES),
        #motif_enrichment
        # expand("motifs/datavis/allmotifs-{condition}-v-{control}-libsizenorm.svg", zip, condition=conditiongroups, control=controlgroups),
        # expand("motifs/datavis/allmotifs-{condition}-v-{control}-spikenorm.svg", zip, condition=conditiongroups_si, control=controlgroups_si)
        expand(expand("motifs/{condition}-v-{control}/{condition}-v-{control}_libsizenorm-{{direction}}-{{category}}-motif_enrichment.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], category=CATEGORIES),
        expand(expand("motifs/{condition}-v-{control}/{condition}-v-{control}_spikenorm-{{direction}}-{{category}}-motif_enrichment.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], category=CATEGORIES),

def plotcorrsamples(wildcards):
    dd = SAMPLES if wildcards.status=="all" else PASSING
    if wildcards.condition=="all":
        if wildcards.norm=="libsizenorm": #condition==all,norm==lib
            return list(dd.keys())
        else: #condition==all,norm==spike
            return [k for k,v in dd.items() if v["spikein"]=="y"]
    elif wildcards.norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in dd.items() if v["group"]==wildcards.control or v["group"]==wildcards.condition]
    else: #condition!=all;norm==spike
        return [k for k,v in dd.items() if (v["group"]==wildcards.control or v["group"]==wildcards.condition) and v["spikein"]=="y"]

rule fastqc_raw:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        "qual_ctrl/fastqc/raw/{sample}"
    threads: config["threads"]
    log: "logs/fastqc/raw/fastqc-raw-{sample}.log"
    shell: """
        (mkdir -p {output}) &> {log}
        (fastqc -o {output} --noextract -t {threads} {input}) &>> {log}
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
        (pigz -f fastq/cleaned/{wildcards.sample}-clean.fastq) &>> {log}
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
        (mkdir -p qual_ctrl/fastqc/cleaned/{wildcards.sample}) &> {log}
        (fastqc -o qual_ctrl/fastqc/cleaned/{wildcards.sample} --noextract -t {threads} {input}) &>> {log}
        """

#align to combined genome with Tophat2, WITHOUT reference transcriptome (i.e., the -G gff)
#(because we don't always have a reference gff and it doesn't make much difference)
rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"]
    output:
        expand(config["tophat2"]["index-path"] + "/" + "{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["tophat2"]["index-path"] + "/" + "{{basename}}.rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2])
    params:
        idx_path = config["tophat2"]["index-path"],
    log: "logs/bowtie2_build.log"
    shell: """
        (bowtie2-build {input.fasta} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand(config["tophat2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".{num}.bt2", num=[1,2,3,4]),
        expand(config["tophat2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2]),
        fastq = "fastq/cleaned/{sample}-clean.fastq.gz"
    output:
        "alignment/{sample}/accepted_hits.bam"
    params:
        idx_path = config["tophat2"]["index-path"],
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
    shell:
        """
        (tophat2 --read-mismatches {params.read_mismatches} --read-gap-length {params.read_gap_length} --read-edit-dist {params.read_edit_dist} -o alignment/{wildcards.sample} --min-anchor-length {params.min_anchor_length} --splice-mismatches {params.splice_mismatches} --min-intron-length {params.min_intron_length} --max-intron-length {params.max_intron_length} --max-insertion-length {params.max_insertion_length} --max-deletion-length {params.max_deletion_length} --num-threads {threads} --max-multihits {params.max_multihits} --library-type fr-firststrand --segment-mismatches {params.segment_mismatches} --no-coverage-search --segment-length {params.segment_length} --min-coverage-intron {params.min_coverage_intron} --max-coverage-intron {params.max_coverage_intron} --min-segment-intron {params.min_segment_intron} --max-segment-intron {params.max_segment_intron} --b2-sensitive {params.idx_path}/{params.basename} {input.fastq}) &> {log}
        """

rule select_unique_mappers:
    input:
        "alignment/{sample}/accepted_hits.bam"
    output:
        temp("alignment/{sample}-unique.bam")
    threads: config["threads"]
    log: "logs/select_unique_mappers/select_unique_mappers-{sample}.log"
    shell: """
        (samtools view -buh -q 50 -@ {threads} {input} | samtools sort -T .{wildcards.sample} -@ {threads} - > {output}) &> {log}
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
    params:
        prefix = lambda wildcards: config["combinedgenome"]["experimental_prefix"] if wildcards.counttype=="counts" else config["combinedgenome"]["spikein_prefix"]
    output:
        plmin = "coverage/{counttype}/{sample}-tss-{counttype}-plmin.bedgraph",
        plus = "coverage/{counttype}/{sample}-tss-{counttype}-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}-tss-{counttype}-minus.bedgraph"
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/get_coverage/get_coverage-{sample}.log"
    shell: """
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.plmin}) &>> {log}
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        """

rule normalize:
    input:
        counts = "coverage/counts/{sample}-tss-counts-{strand}.bedgraph",
        plmin = lambda wildcards: "coverage/counts/" + wildcards.sample + "-tss-counts-plmin.bedgraph" if wildcards.norm=="libsizenorm" else "coverage/sicounts/" + wildcards.sample + "-tss-sicounts-plmin.bedgraph"
    params:
        scalefactor = lambda wildcards: config["spikein-pct"] if wildcards.norm=="spikenorm" else 1
    output:
        normalized = "coverage/{norm}/{sample}-tss-{norm}-{strand}.bedgraph",
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus"
    log: "logs/normalize/normalize-{sample}-{norm}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.plmin} {input.counts} {params.scalefactor} > {output.normalized}) &> {log}
        """

rule get_si_pct:
    input:
        plmin = "coverage/counts/{sample}-tss-counts-plmin.bedgraph",
        SIplmin = "coverage/sicounts/{sample}-tss-sicounts-plmin.bedgraph"
    output:
        temp("qual_ctrl/all/{sample}-spikeincounts.tsv")
    params:
        group = lambda wildcards: SAMPLES[wildcards.sample]["group"]
    log: "logs/get_si_pct/get_si_pct-{sample}.log"
    shell: """
        (echo -e "{wildcards.sample}\t{params.group}\t" $(awk 'BEGIN{{FS=OFS="\t"; ex=0; si=0}}{{if(NR==FNR){{si+=$4}} else{{ex+=$4}}}} END{{print ex+si, ex, si}}' {input.SIplmin} {input.plmin}) > {output}) &> {log}
        """

rule cat_si_pct:
    input:
        expand("qual_ctrl/all/{sample}-spikeincounts.tsv", sample=SAMPLES)
    output:
        "qual_ctrl/all/spikein-counts.tsv"
    log: "logs/cat_si_pct.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_si_pct:
    input:
        "qual_ctrl/all/spikein-counts.tsv"
    output:
        plot = "qual_ctrl/{status}/{status}-spikein-plots.svg",
        stats = "qual_ctrl/{status}/{status}-spikein-stats.tsv"
    params:
        samplelist = lambda wildcards : [k for k,v in SAMPLES.items() if v["spikein"]=="y"] if wildcards.status=="all" else [k for k,v in PASSING.items() if v["spikein"]=="y"],
        conditions = config["comparisons"]["spikenorm"]["conditions"],
        controls = config["comparisons"]["spikenorm"]["controls"],
    script: "scripts/plotsipct.R"

#make 'stranded' genome for datavis purposes
rule make_stranded_genome:
    input:
        exp = config["genome"]["chrsizes"],
        si = config["genome"]["sichrsizes"]
    output:
        exp = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
        si = os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv",
    log: "logs/make_stranded_genome.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.exp} > {output.exp}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.si} > {output.si}) &>> {log}
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
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus} > {output.antisense}) &>> {log}
        """

rule make_stranded_annotations:
    input:
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "{annopath}/stranded/{annotation}-STRANDED.{ext}"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

def selectchrom(wildcards):
    if wildcards.strand in ["plus", "minus"]:
        if wildcards.norm=="sicounts":
            return config["genome"]["sichrsizes"]
        return config["genome"]["chrsizes"]
    if wildcards.norm=="sicounts":
        return os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv"
    return os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"

rule bg_to_bw:
    input:
        bedgraph = "coverage/{norm}/{sample}-tss-{norm}-{strand}.bedgraph",
        chrsizes = selectchrom
    output:
        "coverage/{norm}/{sample}-tss-{norm}-{strand}.bw",
    log : "logs/bg_to_bw/bg_to_bw-{sample}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

rule deeptools_matrix:
    input:
        annotation = lambda wildcards: os.path.dirname(config["annotations"][wildcards.annotation]["path"]) + "/stranded/" + wildcards.annotation + "-STRANDED" + os.path.splitext(config["annotations"][wildcards.annotation]["path"])[1],
        bw = "coverage/{norm}/{sample}-tss-{norm}-{strand}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv"),
        matrix_gz = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv.gz"
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"] + config["annotations"][wildcards.annotation]["binsize"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"] + config["annotations"][wildcards.annotation]["binsize"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}-{norm}-{strand}.log"
    run:
        if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}; pigz -fk {output.matrix}) &> {log}")
        else:
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}; pigz -fk {output.matrix}) &> {log}")

rule melt_matrix:
    input:
        matrix = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv.gz"
    output:
        temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}-melted.tsv.gz")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
    script:
        "scripts/melt_matrix.R"

rule cat_matrices:
    input:
        expand("datavis/{{annotation}}/{{norm}}/{{annotation}}-{sample}-{{norm}}-{{strand}}-melted.tsv.gz", sample=SAMPLES)
    output:
        "datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{annotation}-{norm}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule r_datavis:
    input:
        matrix = "datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{strand}.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{status}_{condition}-v-{control}-{strand}-heatmap-bysample.svg",
        heatmap_group = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{status}_{condition}-v-{control}-{strand}-heatmap-bygroup.svg"
    params:
        samplelist = plotcorrsamples,
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        cluster = lambda wildcards : config["annotations"][wildcards.annotation]["cluster"],
        nclust = lambda wildcards: config["annotations"][wildcards.annotation]["nclusters"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plotHeatmaps.R"

rule union_bedgraph:
    input:
        exp = expand("coverage/{{norm}}/{sample}-tss-{{norm}}-SENSE.bedgraph", sample=SAMPLES)
    output:
        exp = "coverage/{norm}/union-bedgraph-allsamples-{norm}.tsv.gz",
    params:
        names = " ".join(SAMPLES)
    log: "logs/union_bedgraph-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz > {output.exp}) &> {log}
        """

#TODO: plot correlations over multiple binsizes, as well as over annotations
rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-allsamples-{norm}.tsv.gz"
    output:
        "qual_ctrl/{status}/{condition}-v-{control}-tss-{status}-{norm}-correlations.svg"
    params:
        pcount = 0.1,
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule call_tss_peaks:
    input:
        bw = lambda wildcards: "coverage/counts/" + wildcards.sample + "-tss-counts-SENSE.bw" if wildcards.type=="exp" else "coverage/sicounts/" + wildcards.sample + "-tss-sicounts-SENSE.bw"
    output:
        smoothed = expand("peakcalling/{{sample}}-{{type}}-smoothed-bw{bandwidth}-{strand}.bw", strand=["plus","minus"], bandwidth = config["peakcalling"]["bandwidth"]),
        peaks = "peakcalling/{sample}-{type}-allpeaks.narrowPeak"
    params:
        name = lambda wildcards: wildcards.sample + "-" + wildcards.type,
        bandwidth = config["peakcalling"]["bandwidth"],
        window = config["peakcalling"]["local-bg-window"]
    log: "logs/call_tss_peaks/call_tss_peaks-{sample}-{type}.log"
    shell: """
        (python scripts/tss-peakcalling.py -i {input.bw} -n {params.name} -w {params.window} -b {params.bandwidth} -o peakcalling) &> {log}
        """

rule tss_peaks_idr:
    input:
        #NOTE: for now we take the first two samples since the IDR script only takes two
        #change this if we find a better way to aggregate results
        lambda wildcards: ["peakcalling/" + x + "-" + wildcards.type + "-allpeaks.narrowPeak" for x in PASSING if PASSING[x]['group']==wildcards.group][0:2]
    output:
        allpeaks = "peakcalling/{group}-{type}-idrpeaks-all.tsv",
        filtered = "peakcalling/{group}-{type}-idrpeaks-filtered.tsv",
    params:
        idr = int(-125*log2(config["peakcalling"]["idr"]))
    log: "logs/tss_peaks_idr/tss_peaks_idr-{group}-{type}.log"
    shell: """
        idr -s {input} --input-file-type narrowPeak --rank q.value -o {output.allpeaks} -l {log} --plot --peak-merge-method max
        awk -v threshold={params.idr} '$5>threshold || $9=="inf"' peakcalling/{wildcards.group}-{wildcards.type}-idrpeaks-all.tsv > {output.filtered}
        """

rule tss_peaks_to_narrowpeak:
    input:
        "peakcalling/{group}-{type}-idrpeaks-filtered.tsv"
    output:
        "peakcalling/{group}-{type}-idrpeaks.narrowPeak"
    shell: """
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $11, $12, $10}}' {input} | sed -e "s/-minus//g" -e "s/-plus//g" > {output}
        """

rule build_genic_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        windowsize = config["genic-windowsize"]
    log : "logs/build_genic_annotation.log"
    shell: """
        (python scripts/make_genic_annotation.py -t {input.transcripts} -o {input.orfs} -d {params.windowsize} -g {input.chrsizes} -p {output}) &> {log}
        """

rule classify_peaks_genic:
    input:
        peaks = "peakcalling/{group}-exp-idrpeaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/genic/{group}-exp-idrpeaks-genic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | cut --complement -f11-13,15-17 > {output}
        """

rule classify_peaks_intragenic:
    input:
        peaks = "peakcalling/{group}-exp-idrpeaks.narrowPeak",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/intragenic/{group}-exp-idrpeaks-intragenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}} $6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule classify_peaks_antisense:
    input:
        peaks = "peakcalling/{group}-exp-idrpeaks.narrowPeak",
        transcripts = config["genome"]["transcripts"]
    output:
        "peakcalling/antisense/{group}-exp-idrpeaks-antisense.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.transcripts} -wo -S | awk 'BEGIN{{FS=OFS="\t"}} $6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}}$6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule build_convergent_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
    output:
        os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed"
    params:
        max_dist = config["max-convergent-dist"]
    log: "logs/build_convergent_annotation.log"
    shell: """
        (awk -v adist={params.max_dist} 'BEGIN{{FS=OFS="\t"}} $6=="+" {{ if(($3-$2)>adist) print $1, $2, $2+adist, $4, $5, "-" ; else print $0 }} $6=="-" {{if (($3-$2)>adist) print $1, $3-adist, $3, $4, $5, "+"; else print $0}}' {input.transcripts} > {output}) &> {log}
        """

rule classify_peaks_convergent:
    input:
        peaks = "peakcalling/{group}-exp-idrpeaks.narrowPeak",
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/convergent/{group}-exp-idrpeaks-convergent.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}}$6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule build_divergent_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed"
    params:
        max_dist = config["max-divergent-dist"]
    log: "logs/build_divergent_annotation.log"
    shell: """
        (bedtools flank -l {params.max_dist} -r 0 -s -i {input.transcripts} -g {input.chrsizes} | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $3, $4, $5, "-"}} $6=="-"{{print $1, $2, $3, $4, $5, "+"}}' > {output}) &> {log}
        """

rule classify_peaks_divergent:
    input:
        peaks = "peakcalling/{group}-exp-idrpeaks.narrowPeak",
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/divergent/{group}-exp-idrpeaks-divergent.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}}$6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $12-$2+$10}}$6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$13}}' > {output}
        """

rule build_intergenic_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"
    params:
        genic_up = config["genic-windowsize"]
    log: "logs/build_intergenic_annotation.log"
    shell: """
        (bedtools slop -s -l {params.genic_up} -r 0 -i {input.transcripts} -g <(sort -k1,1 {input.chrsizes}) | sort -k1,1 -k2,2n | bedtools complement -i stdin -g <(sort -k1,1 {input.chrsizes}) > {output}) &> {log}
        """

rule classify_peaks_intergenic:
    input:
        peaks = "peakcalling/{group}-exp-idrpeaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/intergenic/{group}-exp-idrpeaks-intergenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.transcripts} {input.orfs} -wa -v | bedtools intersect -a stdin -b {input.genic_anno} -wa -v -s | bedtools intersect -a stdin -b {input.annotation} -wa > {output}
        """

rule peakstats:
    input:
        expand("peakcalling/{category}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
    output:
        table = "peakcalling/{condition}-v-{control}-peaknumbers.tsv",
        size = "peakcalling/{condition}-v-{control}-peaksizes-histogram.svg",
        violin_area = "peakcalling/{condition}-v-{control}-peaksizes-violin-equalarea.svg",
        violin_count = "peakcalling/{condition}-v-{control}-peaksizes-violin-countscaled.svg",
        dist = "peakcalling/{condition}-v-{control}-peakdistances.svg"
    params:
        groups = lambda wildcards: [g for sublist in zip(controlgroups, conditiongroups) for g in sublist] if wildcards.condition=="all" else [wildcards.control, wildcards.condition]
    script:
        "scripts/peakstats.R"

rule combine_tss_peaks:
    input:
        cond = "peakcalling/{condition}-{type}-idrpeaks-filtered.tsv",
        ctrl = "peakcalling/{control}-{type}-idrpeaks-filtered.tsv",
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-peaks.bed"
    shell: """
        sort -k1,1 -k2,2n {input.cond} | bedtools multiinter -i stdin <(sort -k1,1 -k2,2n {input.ctrl}) -cluster | cut -f1-3 | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule map_counts_to_peaks:
    input:
        bed = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-peaks.bed",
        bg = lambda wildcards: "coverage/counts/" + wildcards.sample + "-tss-counts-SENSE.bedgraph" if wildcards.type=="exp" else "coverage/sicounts/" + wildcards.sample + "-tss-sicounts-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/{sample}-{type}-allpeakcounts.tsv")
    log: "logs/map_counts_to_peaks/map_counts_to_peaks-{condition}-v-{control}-{sample}-{type}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-"$2"-"$3, $4}}' &> {output}) &> {log}
        """

def getsamples(ctrl, cond):
    return [k for k,v in PASSING.items() if (v["group"]==ctrl or v["group"]==cond)]

rule get_peak_counts:
    input:
        lambda wildcards : ["diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/" + x + "-" + wildcards.type + "-allpeakcounts.tsv" for x in getsamples(wildcards.control, wildcards.condition)]
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-peak-counts.tsv"
    params:
        n = lambda wildcards: 2*len(getsamples(wildcards.control, wildcards.condition)),
        names = lambda wildcards: "\t".join(getsamples(wildcards.control, wildcards.condition))
    log: "logs/get_peak_counts/get_peak_counts-{condition}-v-{control}-{type}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule call_de_peaks:
    input:
        expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-peak-counts.tsv",
        sicounts = lambda wildcards: "diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-si-peak-counts.tsv" if wildcards.norm=="spikenorm" else "diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-exp-peak-counts.tsv"
    params:
        samples = lambda wildcards : getsamples(wildcards.control, wildcards.condition),
        groups = lambda wildcards : [PASSING[x]["group"] for x in getsamples(wildcards.control, wildcards.condition)],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    output:
        results = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv",
        #need to write out norm counts here or just in the total qc?
        normcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-counts-sfnorm-{norm}.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-counts-rlog-{norm}.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-qcplots-{norm}.svg"
    script:
        "scripts/call_de_peaks.R"

rule separate_de_peaks:
    input:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv",
    output:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-up.tsv",
        unchanged = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-unchanged.tsv",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-down.tsv",
    params:
        fdr = -log10(config["deseq"]["fdr"])
    log: "logs/separate_de_peaks/separate_de_peaks-{condition}-v-{control}-{norm}.log"
    shell: """
        (awk -v afdr={params.fdr} 'BEGIN{{FS=OFS="\t"}} NR==1{{print > "{output.up}"; print > "{output.unchanged}"; print > "{output.down}" }} NR>1 && $10>afdr && $7>0 {{print > "{output.up}"}} NR>1 && $10>afdr && $7<0 {{print > "{output.down}"}} NR>1 && $10<=afdr {{print > "{output.unchanged}"}}' {input}) &> {log}
        """

rule de_peaks_to_bed:
    input:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.tsv",
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
    log: "logs/de_peaks_to_bed/de_peaks_to_bed-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (tail -n +2 {input} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, $7":"$11, $3}}' > {output}) &> {log}
        """

rule get_de_intragenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.bed",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-results-{norm}-all-intragenic.tsv"
    log: "logs/get_putative_intragenic/get_putative_intragenic-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}} $6=="-"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\tORF_start\tORF_end\tORF_name\tdist_ATG_to_peak") - > {output}) &> {log}
        """

rule get_de_intragenic_frequency:
    input:
        orfs = config["genome"]["orf-annotation"],
        intrabed = "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic.bed"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
    shell: """
        bedtools intersect -a {input.orfs} -b {input.intrabed} -c -s > {output}
        """

rule plot_de_intragenic_frequency:
    input:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-{norm}-{direction}-freqperORF.svg"
    log: "logs/plot_intragenic_frequency/plot_intragenic_frequency-{condition}-v-{control}-{norm}-{direction}.log"
    script: "scripts/intrafreq.R"

rule get_de_antisense:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.bed",
        transcripts = config["genome"]["transcripts"],
        totalresults = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv"
    output:
        "diff_exp/{condition}-v-{control}/antisense/{condition}-v-{control}-results-{norm}-all-antisense.tsv"
    log : "logs/get_de_antisense/get_de_antisense-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcripts} -wo -S | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}} $6=="-"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name\tdist_peak_to_senseTSS") - > {output}) &> {log}
        """

rule get_de_genic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv"
    output:
        "diff_exp/{condition}-v-{control}/genic/{condition}-v-{control}-results-{norm}-all-genic.tsv"
    log : "logs/get_de_genic/get_de_genic-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10}} $6=="-"{{print $4, $8, $9, $10}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name") - > {output}) &> {log}
        """

rule get_de_intergenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv"
    output:
        "diff_exp/{condition}-v-{control}/intergenic/{condition}-v-{control}-results-{norm}-all-intergenic.tsv"
    log : "logs/get_de_intergenic/get_de_intergenic-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcripts} {input.orfs} -wa -v | bedtools intersect -a stdin -b {input.genic_anno} -wa -v -s | bedtools intersect -a stdin -b {input.annotation} -wo | awk 'BEGIN{{FS=OFS="\t"}}$6=="+"{{print $4, $8, $9}}$6=="-"{{print $4, $8, $9}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\tregion_start\tregion_end") - > {output}) &> {log}
        """

rule get_intra_orfs:
    input:
        peaks = "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-de-clusters-{norm}-{direction}-intragenic.tsv",
        fasta = config["genome"]["fasta"]
    output:
        "diff_exp/{condition}-v-{control}/intragenic/intragenic-orfs/{condition}-v-{control}-{norm}-{direction}-intragenic-orfs.tsv"
    params:
        max_upstr_atgs = config["max-upstr-atgs"],
        max_search_dist = 2000
    log: "logs/get_intra_orfs/get_intra_orfs-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (python scripts/find_intra_orfs.py -p {input.peaks} -f {input.fasta} -m {params.max_search_dist} -a {params.max_upstr_atgs} -o {output}) &> {log}
        """

rule get_de_convergent:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.bed",
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv"
    output:
        "diff_exp/{condition}-v-{control}/convergent/{condition}-v-{control}-results-{norm}-all-convergent.tsv"
    log : "logs/get_de_convergent/get_de_convergent-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}} $6=="-"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name\tdist_peak_to_senseTSS") - > {output}) &> {log}
        """

rule get_de_divergent:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.bed",
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv"
    output:
        "diff_exp/{condition}-v-{control}/divergent/{condition}-v-{control}-results-{norm}-all-divergent.tsv"
    log : "logs/get_de_divergent/get_de_divergent-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}} $6=="-"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name\tdist_peak_to_senseTSS") - > {output}) &> {log}
        """

rule separate_sig_de:
    input:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-all-{category}.tsv"
    output:
        up = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-up-{category}.tsv",
        unchanged = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-unchanged-{category}.tsv",
        down = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-down-{category}.tsv"
    params:
        fdr = -log10(config["deseq"]["fdr"])
    log: "logs/separate_sig_de/separate_sig_de-{condition}-v-{control}-{norm}-{category}.log"
    shell: """
        awk -v afdr={params.fdr} 'BEGIN{{FS=OFS="\t"}}NR==1{{print > "{output.up}"; print > "{output.unchanged}"; print > "{output.down}"}} NR>1 && $7>0 && $8>afdr {{print > "{output.up}"}} NR>1 && $7<0 && $8>afdr {{print > "{output.down}"}} NR>1 && $8<=afdr{{print > "{output.unchanged}"}}' {input}
        """

rule get_de_category_bed:
    input:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.tsv"
    output:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed"
    log: "logs/get_category_bed/get_category_bed-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (tail -n +2 {input} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, $7":"$8, $3}}' | sort -k1,1 -k2,2n  > {output}) &> {log}
        """

#TODO: account for double-counted peaks when a peak overlaps more than one annotation (more than one genic region, for example)
rule summarise_de_results:
    input:
        total = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv",
        genic = "diff_exp/{condition}-v-{control}/genic/{condition}-v-{control}-results-{norm}-all-genic.tsv",
        intragenic = "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-results-{norm}-all-intragenic.tsv",
        antisense = "diff_exp/{condition}-v-{control}/antisense/{condition}-v-{control}-results-{norm}-all-antisense.tsv",
        convergent = "diff_exp/{condition}-v-{control}/convergent/{condition}-v-{control}-results-{norm}-all-convergent.tsv",
        divergent = "diff_exp/{condition}-v-{control}/divergent/{condition}-v-{control}-results-{norm}-all-divergent.tsv",
        intergenic = "diff_exp/{condition}-v-{control}/intergenic/{condition}-v-{control}-results-{norm}-all-intergenic.tsv",
    params:
        lfc = log2(config["deseq"]["fold-change-threshold"]),
        alpha = config["deseq"]["fdr"]
    output:
        summary = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{norm}-diffexp-summary.svg",
        maplot = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{norm}-diffexp-maplot.svg",
        volcano = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{norm}-diffexp-volcano.svg",
    script: "scripts/de_summary.R"

rule genic_v_class:
    input:
        genic = "diff_exp/{condition}-v-{control}/genic/{condition}-v-{control}-results-{norm}-all-genic.tsv",
        intragenic = "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-results-{norm}-all-intragenic.tsv",
        antisense = "diff_exp/{condition}-v-{control}/antisense/{condition}-v-{control}-results-{norm}-all-antisense.tsv",
        convergent = "diff_exp/{condition}-v-{control}/convergent/{condition}-v-{control}-results-{norm}-all-convergent.tsv",
        divergent = "diff_exp/{condition}-v-{control}/divergent/{condition}-v-{control}-results-{norm}-all-divergent.tsv",
    params:
        path = "diff_exp/{condition}-v-{control}/genic_v_class/"
    output:
        figure = "diff_exp/{condition}-v-{control}/genic_v_class/{condition}-v-{control}-{norm}-genic-v-class.svg",
        tables = expand("diff_exp/{{condition}}-v-{{control}}/genic_v_class/{{condition}}-v-{{control}}-{{norm}}-genic-v-{ttype}.tsv", ttype=["intragenic", "antisense", "convergent", "divergent"])
    script: "scripts/classvgenic.R"

#get all motif names from motif databases, cleaning nasty characters in some motif names
MOTIFS = subprocess.run(args="meme2meme " + " ".join(config["motifs"]["databases"]) + " | grep -e '^MOTIF' | cut -d ' ' -f2 | sed 's/\//_/g; s/&/_/g; s/{/[/g; s/}/]/g' ", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()

#run fimo in parallel for each motif for speed
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
        fimo --motif {wildcards.motif} --bgfile <(fasta-get-markov {input.fasta}) --parse-genomic-coord --thresh {params.alpha} --text <(meme2meme {input.motif_db}) {input.fasta} | sed -e 's/\//_/g; s/&/_/g; s/{{/[/g; s/}}/]/g' | tee {output.tsv} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $3, $4-1, $5, $1, -log($8)/log(10), $6, $2, $10}}' | sort -k1,1 -k2,2n > {output.bed}
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
rule get_upstream_motifs:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed",
        chrsizes = config["genome"]["chrsizes"],
        motifs = "motifs/allmotifs.bed"
    output:
        "motifs/{condition}-v-{control}/{condition}-v-{control}_{norm}-{direction}-{category}-motifs.tsv.gz"
    params:
        upstr = config["motifs"]["upstream"],
        dnstr = config["motifs"]["enrichment-downstream"]
    log: "logs/get_upstream_motifs/get_upstream_motifs-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (uniq {input.peaks} | bedtools flank -l {params.upstr} -r 0 -s -i stdin -g {input.chrsizes} | bedtools slop -l 0 -r {params.dnstr} -s -i stdin -g {input.chrsizes} | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b {input.motifs} -sorted -wao | awk 'BEGIN{{FS="\t|:"; OFS="\t"}}{{print $1, $4, $5, $6, $11, $14, $9, $10, $12}}' | cat <(echo -e "chrom\ttss_peak_id\tpeak_lfc\tpeak_logpadj\tmotif_id\tmotif_alt_id\tmotif_start\tmotif_end\tmotif_logpadj") - | pigz -f > {output}) &> {log}
        """

rule test_motif_enrichment:
    input:
        fimo_pos = "motifs/{condition}-v-{control}/{condition}-v-{control}_{norm}-{direction}-{category}-motifs.tsv.gz",
        fimo_neg = "motifs/{condition}-v-{control}/{condition}-v-{control}_{norm}-unchanged-{category}-motifs.tsv.gz"
    params:
        pval_cutoff = config["motifs"]["fimo-pval"],
        alpha= config["motifs"]["enrichment-fdr"],
        direction = lambda wildcards: "upregulated" if wildcards.direction=="up" else "downregulated"
    output:
        tsv = "motifs/{condition}-v-{control}/{condition}-v-{control}_{norm}-{direction}-{category}-motif_enrichment.tsv",
        plot = "motifs/{condition}-v-{control}/{condition}-v-{control}_{norm}-{direction}-{category}-motif_enrichment.svg",
    script: "scripts/motif_enrichment.R"

rule get_motif_coverage:
    input:
        bed = "motifs/.{motif}.bed", #this is sorted when created
        chrsizes = config["genome"]["chrsizes"]
    output:
        bg = "motifs/coverage/{motif}.bedgraph",
        bw = "motifs/coverage/{motif}.bw",
    shell: """
        cut -f1-6 {input.bed} | bedtools genomecov -bga -i stdin -g {input.chrsizes} | sort -k1,1 -k2,2n > {output.bg}
        bedGraphToBigWig {output.bg} {input.chrsizes} {output.bw}
        """

rule motif_matrix:
    input:
        annotation = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed",
        bw = "motifs/coverage/{motif}.bw"
    output:
        dtfile = temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.mat"),
        matrix = temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv"),
        matrix_gz = "motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv.gz",
    params:
        refpoint = "TSS",
        upstream = config["motifs"]["upstream"] + config["motifs"]["binsize"],
        dnstream = config["motifs"]["freq-downstream"] + config["motifs"]["binsize"],
        binsize = config["motifs"]["binsize"],
        sort = "keep",
        binstat = "sum"
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{motif}-{condition}-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --averageTypeBins {params.binstat} -p {threads}) &> {log}
        pigz -fk {output.matrix}
        """

rule melt_motif_matrix:
    input:
        matrix = "motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv.gz",
    output:
        temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}-melted.tsv.gz"),
    params:
        refpoint = "TSS",
        binsize = config["motifs"]["binsize"],
        upstream = config["motifs"]["upstream"],
    script:
        "scripts/melt_motif_matrix.R"

rule cat_motif_matrices:
    input:
        expand("motifs/datavis/{motif}_{{condition}}-v-{{control}}_{{norm}}-{direction}-peaks-{category}-melted.tsv.gz", category=CATEGORIES, motif=MOTIFS, direction=["up","down","unchanged"]),
    output:
        "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{condition}-{control}-{norm}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_motif_freq:
    input:
        "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.tsv.gz"
    output:
        "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.svg"
    script: "scripts/motif_metagenes.R"

#for peaks are double-counted; only keep one sequence if two are overlapping
# rule get_peak_sequences_nooverlap:
#     input:
#         peaks = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed",
#         chrsizes = config["genome"]["chrsizes"],
#         fasta = config["genome"]["fasta"]
#     output:
#         "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}-nooverlap.fa"
#     params:
#         upstr = config["meme-chip"]["upstream-dist"],
#         dnstr = config["meme-chip"]["downstream-dist"]
#     log: "logs/get_peak_sequences/get_peak_sequences-{condition}-v-{control}-{norm}-{direction}-{category}.log"
#     shell: """
#         (uniq {input.peaks} | bedtools flank -l {params.upstr} -r 0 -s -i stdin -g {input.chrsizes} | bedtools slop -l 0 -r {params.dnstr} -s -i stdin -g {input.chrsizes} | awk 'BEGIN{{FS=OFS="\t"}}{{$6=="-" ? $1=$1"-minus":$1=$1"-plus"}}{{print $0}}' | sort -k1,1 -k2,2n | bedtools spacing -i stdin | awk 'BEGIN{{FS=OFS="\t"}} $7!=0{{print $1, $2, $3, $4, $5, $6}}' | sed -e 's/-minus//g;s/-plus//g' | bedtools getfasta -name -s -fi {input.fasta} -bed stdin > {output}) &> {log}
#         """
# rule fimo:
#     input:
#         fa = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}-fimo.fa",
#         dbs = config["meme-chip"]["motif-databases"]
#     params:
#         pval = config["meme-chip"]["fimo-pval"],
#         qval = config["meme-chip"]["fimo-qval"],
#     output:
#         txt = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-{norm}-{direction}-{category}-fimo/fimo.txt",
#         gff_all = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-{norm}-{direction}-{category}-fimo/{condition}-v-{control}-{norm}-{direction}-{category}-fimo_all.gff",
#         gff_filtered = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-{norm}-{direction}-{category}-fimo/{condition}-v-{control}-{norm}-{direction}-{category}-fimo_filtered.gff",
#     shell: """
#         fimo --bgfile <(fasta-get-markov {input.fa}) --oc diff_exp/{wildcards.condition}-v-{wildcards.control}/{wildcards.category}/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-{wildcards.direction}-{wildcards.category}-fimo --thresh {params.pval} --parse-genomic-coord <(meme2meme {input.dbs}) {input.fa}
#         sed 's/peak_[[:digit:]]\+_//g' diff_exp/{wildcards.condition}-v-{wildcards.control}/{wildcards.category}/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-{wildcards.direction}-{wildcards.category}-fimo/fimo.gff | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{$4=$4+1; $5=$5+1}}{{print $0}}' | awk -v alpha={params.qval} 'BEGIN{{FS="qvalue= |;"}} {{print $0 > "{output.gff_all}"}} NR==1 || $6<alpha {{print $0 > "{output.gff_filtered}"}}'
#         """

# rule meme_chip:
#     input:
#         seq = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.fa",
#         background = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-unchanged-{category}.fa",
#         dbs = config["meme-chip"]["motif-databases"]
#     output:
#         "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-{norm}-{direction}-{category}-motifs/meme-chip.html"
#     params:
#         dbs = ["-db " + x for x in config["meme-chip"]["motif-databases"]],
#         ccut = config["meme-chip"]["max-frag-size"],
#         mode = config["meme-chip"]["meme-mode"],
#         nmotifs = config["meme-chip"]["meme-nmotifs"],
#     conda: "envs/meme_chip.yaml"
#     shell: """
#         meme-chip {input.seq} -neg {input.background} -oc diff_exp/{wildcards.condition}-v-{wildcards.control}/{wildcards.category}/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-{wildcards.direction}-{wildcards.category}-motifs {params.dbs} -ccut {params.ccut} -meme-mod {params.mode} -meme-nmotifs {params.nmotifs} -nmeme 10000 -meme-maxsize 600000
#         """
