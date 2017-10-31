#!/usr/bin/env python
import os
from math import log2

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

CATEGORIES = ["genic", "intragenic", "intergenic", "antisense", "convergent", "divergent"]

localrules: all,
    bowtie2_build,
    get_si_pct,
    cat_si_pct,
    plot_si_pct,
    make_stranded_annotations,
    map_counts_to_peaks,
    get_peak_counts,
    de_clusters_to_bed,
    extract_base_cluster_dist,
    map_counts_to_genic,
    get_genic_counts,
    get_putative_intragenic,
    get_intragenic_frequency,
    # plot_intragenic_frequency,
    get_putative_antisense,
    build_genic_annotation,
    get_putative_genic,
    build_intergenic_annotation,
    get_putative_intergenic,
    # get_intra_orfs,
    build_convergent_annotation,
    get_putative_convergent,
    build_divergent_annotation,
    get_putative_divergent,
    get_category_bed,
    # get_peak_sequences,
    # meme_chip,
    # class_v_genic,
rule all:
    input:
        #FastQC
        expand("qual_ctrl/fastqc/raw/{sample}", sample=SAMPLES),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/bw/{sample}-tss-{norm}-{strand}.bw", norm=["spikenorm","libsizenorm", "counts", "sicounts"], sample=SAMPLES, strand=["SENSE","ANTISENSE","plus","minus"]),
        #datavis
        # expand(expand("datavis/{{annotation}}/libsizenorm/tss-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}-{{strand}}-heatmap-bygroup.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], status=["all","passing"], strand=["SENSE","ANTISENSE"]),
        # expand(expand("datavis/{{annotation}}/spikenorm/tss-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}-{{strand}}-heatmap-bygroup.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), annotation=config["annotations"], status=["all","passing"], strand=["SENSE","ANTISENSE"]),
        #quality control
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all", "passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-tss-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status = ["all", "passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-tss-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status = ["all", "passing"]),
        #find intragenic ORFs
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intragenic-orfs/{condition}-v-{control}-libsizenorm-{{direction}}-intragenic-orfs.tsv", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"]),
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intragenic-orfs/{condition}-v-{control}-spikenorm-{{direction}}-intragenic-orfs.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"]),
        #MEME-ChIP
        # expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-spikenorm-{{direction}}-{{category}}-motifs/index.html", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"], category = CATEGORIES),
        # expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-libsizenorm-{{direction}}-{{category}}-motifs/index.html", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"], category = CATEGORIES),
        # expand(expand("diff_exp/{condition}-v-{control}/{{type}}/{{type}}-v-genic/{condition}-v-{control}-{{type}}-v-genic-spikenorm.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), type=["antisense", "convergent", "divergent", "intragenic"]),
        # expand(expand("diff_exp/{condition}-v-{control}/{{type}}/{{type}}-v-genic/{condition}-v-{control}-{{type}}-v-genic-libsizenorm.tsv", zip, condition=conditiongroups, control=controlgroups), type=["antisense", "convergent", "divergent", "intragenic"]),
        #peakcalling on all samples
        expand("peakcalling/{sample}-exp-allpeaks.narrowPeak", sample=SAMPLES),
        expand("peakcalling/{sample}-si-allpeaks.narrowPeak", sample=sisamples),
        #IDR for all groups which have at least two passing samples
        expand("peakcalling/{group}-exp-idrpeaks.{fmt}", group = validgroups, fmt=["tsv", "narrowPeak"]),
        expand("peakcalling/{group}-si-idrpeaks.{fmt}", group = validgroups_si, fmt=["tsv", "narrowPeak"]),
        #classify peaks into categories
        expand("peakcalling/{category}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
        "peakcalling/peakdistances.svg",
        #combine called peaks for conditions vs control
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-peaks.bed", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-si-peaks.bed", zip, condition=conditiongroups_si, control=controlgroups_si),
        #differential expression of peaks
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-peak-counts.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-si-peak-counts.tsv", zip, condition=conditiongroups_si, control=controlgroups_si),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-libsizenorm-all.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-spikenorm-all.tsv", zip, condition=conditiongroups_si, control=controlgroups_si),
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-libsizenorm-{{dir}}.{{fmt}}", zip, condition=conditiongroups, control=controlgroups), dir=["up","down"], fmt=["tsv","bed"]),
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-spikenorm-{{dir}}.{{fmt}}", zip, condition=conditiongroups_si, control=controlgroups_si), dir=["up","down"], fmt=["tsv","bed"]),
        #categorize DE peaks
        expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-results-libsizenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups, control=controlgroups), direction = ["up","down"], category=CATEGORIES),
        expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-results-spikenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up","down"], category=CATEGORIES),
        #intragenic frequency per ORF
        expand(expand("diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-libsizenorm-{{direction}}-freqperORF.svg", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"]),
        expand(expand("diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-spikenorm-{{direction}}-freqperORF.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"]),

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
        expand("../genome/bowtie2_indexes/{basename}.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2,3,4]),
        expand("../genome/bowtie2_indexes/{basename}.rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2])
    params:
        name = config["combinedgenome"]["name"]
    log: "logs/bowtie2_build.log"
    shell: """
        (bowtie2-build {input.fasta} ../genome/bowtie2_indexes/{params.name}) &> {log}
        """

rule align:
    input:
        expand("../genome/bowtie2_indexes/{basename}.{num}.bt2", basename=config["combinedgenome"]["name"], num = [1,2,3,4]),
        expand("../genome/bowtie2_indexes/{basename}.rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2]),
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
    shell:
        """
        (tophat2 --read-mismatches {params.read_mismatches} --read-gap-length {params.read_gap_length} --read-edit-dist {params.read_edit_dist} -o alignment/{wildcards.sample} --min-anchor-length {params.min_anchor_length} --splice-mismatches {params.splice_mismatches} --min-intron-length {params.min_intron_length} --max-intron-length {params.max_intron_length} --max-insertion-length {params.max_insertion_length} --max-deletion-length {params.max_deletion_length} --num-threads {threads} --max-multihits {params.max_multihits} --library-type fr-firststrand --segment-mismatches {params.segment_mismatches} --no-coverage-search --segment-length {params.segment_length} --min-coverage-intron {params.min_coverage_intron} --max-coverage-intron {params.max_coverage_intron} --min-segment-intron {params.min_segment_intron} --max-segment-intron {params.max_segment_intron} --b2-sensitive ../genome/bowtie2_indexes/{params.basename} {input.fastq}) &> {log}
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
        SIplmin = "coverage/sicounts/{sample}-tss-sicounts-plmin.bedgraph",
        SIpl = "coverage/sicounts/{sample}-tss-sicounts-plus.bedgraph",
        SImin = "coverage/sicounts/{sample}-tss-sicounts-minus.bedgraph",
        plmin = "coverage/counts/{sample}-tss-counts-plmin.bedgraph",
        plus = "coverage/counts/{sample}-tss-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-tss-counts-minus.bedgraph"
    params:
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
        si_prefix = config["combinedgenome"]["spikein_prefix"]
    log: "logs/get_coverage/get_coverage-{sample}.log"
    shell: """
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SIplmin}) &> {log}
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SIpl}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SImin}) &>> {log}
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.plmin}) &>> {log}
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        """

rule normalize:
    input:
        plus = "coverage/counts/{sample}-tss-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-tss-counts-minus.bedgraph",
        plmin = "coverage/counts/{sample}-tss-counts-plmin.bedgraph",
        SIplmin = "coverage/sicounts/{sample}-tss-sicounts-plmin.bedgraph"
    output:
        spikePlus = "coverage/spikenorm/{sample}-tss-spikenorm-plus.bedgraph",
        spikeMinus = "coverage/spikenorm/{sample}-tss-spikenorm-minus.bedgraph",
        libnormPlus = "coverage/libsizenorm/{sample}-tss-libsizenorm-plus.bedgraph",
        libnormMinus = "coverage/libsizenorm/{sample}-tss-libsizenorm-minus.bedgraph"
    params:
        scalefactor = config["spikein-pct"]
    log: "logs/normalize/normalize-{sample}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.SIplmin} {input.plus} {params.scalefactor} > {output.spikePlus}) &> {log}
        (bash scripts/libsizenorm.sh {input.SIplmin} {input.minus} {params.scalefactor} > {output.spikeMinus}) &>> {log}
        (bash scripts/libsizenorm.sh {input.plmin} {input.plus} 1 > {output.libnormPlus}) &>> {log}
        (bash scripts/libsizenorm.sh {input.plmin} {input.minus} 1 > {output.libnormMinus}) &>> {log}
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
        "../genome/annotations/stranded/{annotation}-STRANDED.bed"
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
        "coverage/{norm}/bw/{sample}-tss-{norm}-{strand}.bw",
    log : "logs/bg_to_bw/bg_to_bw-{sample}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
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
    run:
        if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
        else:
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")

rule gzip_deeptools_matrix:
    input:
        matrix = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv"
    output:
        "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv.gz"
    shell: """
        pigz -f {input}
        """

rule melt_matrix:
    input:
        matrix = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv.gz"
    output:
        temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}-melted.tsv.gz")
    params:
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"]
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
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
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

rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-allsamples-{norm}.tsv.gz"
    output:
        "qual_ctrl/{status}/{condition}-v-{control}-tss-{norm}-correlations.svg"
    params:
        pcount = 0.1,
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule call_tss_peaks:
    input:
        bw = lambda wildcards: "coverage/counts/bw/" + wildcards.sample + "-tss-counts-SENSE.bw" if wildcards.type=="exp" else "coverage/sicounts/bw/" + wildcards.sample + "-tss-sicounts-SENSE.bw"
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
        "peakcalling/{group}-{type}-idrpeaks.tsv"
    params:
        idr = config["peakcalling"]["idr"]
    log: "logs/tss_peaks_idr/tss_peaks_idr-{group}-{type}.log"
    shell: """
        idr -s {input} --input-file-type narrowPeak --rank q.value -o {output} -l {log} -i {params.idr} --plot --peak-merge-method max
        """

rule tss_peaks_to_narrowpeak:
    input:
        "peakcalling/{group}-{type}-idrpeaks.tsv"
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
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"
    output:
        "peakcalling/intergenic/{group}-exp-idrpeaks-intergenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.annotation} -wa > {output}
        """

rule peakstats:
    input:
        expand("peakcalling/{category}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
    output:
        table = "peakcalling/peaknumbers.tsv",
        size = "peakcalling/peaksizes.svg",
        dist = "peakcalling/peakdistances.svg"
    params:
        groups = [g for sublist in zip(controlgroups, conditiongroups) for g in sublist]
    script:
        "scripts/peakstats.R"

rule combine_tss_peaks:
    input:
        cond = "peakcalling/{condition}-{type}-idrpeaks.tsv",
        ctrl = "peakcalling/{control}-{type}-idrpeaks.tsv",
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
        alpha = config["deseq"]["fdr"]
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
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-down.tsv",
    params:
        fdr = config["deseq"]["fdr"]
    shell: """
        awk -v afdr={params.fdr} 'BEGIN{{FS=OFS="\t"}} NR==1{{print > "{output.up}"; print > "{output.down}" }} NR>1 && $7<afdr && $3>0 {{print > "{output.up}"}} NR>1 && $7<afdr && $3<0 {{print > "{output.down}"}}' {input}
        """

#NOTE: column 6 for down tables vs column 5 for up tables is to the negative sign in the fold-change
rule de_peaks_to_bed:
    input:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-up.tsv",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-down.tsv",
    output:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-up.bed",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-down.bed",
    log: "logs/de_peaks_to_bed/de_peaks_to_bed-{condition}-v-{control}-{norm}.log"
    shell: """
        (tail -n +2 {input.up} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $3":"(-log($7)/log(10))}}' | awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1, $3, $4, "up_"NR, $5, "+"}} $2=="minus"{{print $1, $3, $4, "up_"NR, $5, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.up}) &> {log}
        (tail -n +2 {input.down} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $3":"(-log($7)/log(10))}}' | awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1, $3, $4, "down_"NR, $6, "+"}} $2=="minus"{{print $1, $3, $4, "down_"NR, $6, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.down}) &>> {log}
        """

#TODO: add a cat statement to add a header for all 'class' tsv (make sure to check class to bed afterwards)
rule get_de_intragenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic.tsv"
    log: "logs/get_putative_intragenic/get_putative_intragenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}}' | sort -k10,10nr | cat <(echo -e "chrom\tstrand\tpeak_start\tpeak_end\tpeak_name\tORF_start\tORF_end\tORF_name\tpeak_lfc\tpeak_significance\tdist_atg_to_peak") - > {output}) &> {log}
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
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        transcripts = config["genome"]["transcripts"]
    output:
        "diff_exp/{condition}-v-{control}/antisense/{condition}-v-{control}-results-{norm}-{direction}-antisense.tsv"
    log : "logs/get_putative_antisense/get_putative_antisense-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcripts} -wo -S | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}}' | sort -k10,10nr | cat <(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance\tdist_peak_to_senseTSS") - > {output}) &> {log}
        """


rule get_de_genic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/genic/{condition}-v-{control}-results-{norm}-{direction}-genic.tsv"
    log : "logs/get_putative_genic/get_putative_genic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6}}' | sort -k10,10nr | cat <(echo -e "chrom\tstrand\tpeak_start\tpeak_end\tpeak_name\tgenic_start\tgenic_end\tgenic_name\tpeak_lfc\tpeak_significance") - > {output}) &> {log}
        """

rule get_de_intergenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/intergenic/{condition}-v-{control}-results-{norm}-{direction}-intergenic.tsv"
    log : "logs/get_putative_intergenic/get_putative_intergenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo | awk 'BEGIN{{FS="\t|:";OFS="\t"}}{{print $1, $7, $2, $3, $4, $9, $10, ".", $5, $6}}'| sort -k10,10nr | cat <(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\tregion_start\tregion_end\tregion_name\tpeak_lfc\tpeak_significance") - > {output}) &> {log}
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
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/convergent/{condition}-v-{control}-results-{norm}-{direction}-convergent.tsv"
    log : "logs/get_putative_convergent/get_putative_convergent-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}}' | sort -k10,10nr | cat <(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance\tdist_peak_to_senseTSS") - > {output}) &> {log}
        """


rule get_de_divergent:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/divergent/{condition}-v-{control}-results-{norm}-{direction}-divergent.tsv"
    log : "logs/get_putative_divergent/get_putative_divergent-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}}' | sort -k10,10nr | cat <(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance\tdist_peak_to_senseTSS") - > {output}) &> {log}
        """

rule get_de_category_bed:
    input:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.tsv"
    output:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed"
    log: "logs/get_category_bed/get_category_bed-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (tail -n +2 {input} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $3, $4, $5, $10, $2}}' | sort -k1,1 -k2,2n  > {output}) &> {log}
        """

rule get_peak_sequences:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-cluster-results-{norm}-{direction}-{category}.bed",
        chrsizes = config["genome"]["chrsizes"],
        fasta = config["genome"]["fasta"]
    output:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-cluster-results-{norm}-{direction}-{category}.fa"
    params:
        upstr = config["meme-chip"]["upstream-dist"],
        dnstr = config["meme-chip"]["downstream-dist"]
    log: "logs/get_peak_sequences/get_peak_sequences-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (bedtools slop -l {params.upstr} -r {params.dnstr} -s -i {input.peaks} -g {input.chrsizes} | bedtools getfasta -name -s -fi {input.fasta} -bed stdin > {output}) &> {log}
        """

rule meme_chip:
    input:
        seq = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-cluster-results-{norm}-{direction}-{category}.fa",
        db = config["meme-chip"]["motif-database"]
    output:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-{norm}-{direction}-{category}-motifs/index.html"
    params:
        ccut = config["meme-chip"]["max-frag-size"],
        mode = config["meme-chip"]["meme-mode"],
        nmotifs = config["meme-chip"]["meme-nmotifs"],
    #threads: config["threads"]
    shell: """
        meme-chip {input.seq} -oc diff_exp/{wildcards.condition}-v-{wildcards.control}/{wildcards.category}/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-{wildcards.direction}-{wildcards.category}-motifs -db {input.db} -ccut {params.ccut} -meme-mod {params.mode} -meme-nmotifs {params.nmotifs} -meme-p 2
        """

rule class_v_genic:
    input:
        pclass_up = "diff_exp/{condition}-v-{control}/{type}/{condition}-v-{control}-cluster-results-{norm}-up-{category}.tsv",
        pclass_dn = "diff_exp/{condition}-v-{control}/{type}/{condition}-v-{control}-cluster-results-{norm}-down-{category}.tsv",
        genic = "diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-results-{norm}-all.tsv",
    output:
        scatter_text = "diff_exp/{condition}-v-{control}/{type}/{type}-v-genic/{condition}-v-{control}-{type}-v-genic-{norm}-scattertext.svg" ,
        scatter_nolabel = "diff_exp/{condition}-v-{control}/{type}/{type}-v-genic/{condition}-v-{control}-{type}-v-genic-{norm}-scatternotext.svg",
        table = "diff_exp/{condition}-v-{control}/{type}/{type}-v-genic/{condition}-v-{control}-{type}-v-genic-{norm}.tsv"
    script:
        "scripts/class_v_genic.R"
