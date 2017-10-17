#!/usr/bin/env python
import os
from math import log2

configfile: "config.yaml"

SAMPLES = config["samples"]
name_string = " ".join(SAMPLES)
PASSING = {k:v for (k,v) in SAMPLES.items() if v["pass-qc"] == "pass"}
pass_string = " ".join(PASSING)

controlgroups = config["comparisons"]["libsizenorm"]["controls"]
conditiongroups = config["comparisons"]["libsizenorm"]["conditions"]
controlgroups_si = config["comparisons"]["spikenorm"]["controls"]
conditiongroups_si = config["comparisons"]["spikenorm"]["conditions"]

CATEGORIES = ["genic", "intragenic", "intergenic", "antisense", "convergent", "divergent"]

localrules: all,
    bowtie2_build,
    get_si_pct,
    cat_si_pct,
    plot_si_pct,
    make_stranded_annotations,
    separate_de_results,
    de_bases_to_bed,
    merge_de_bases_to_clusters,
    map_counts_to_clusters,
    get_cluster_counts,
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
        # expand("coverage/{norm}/bw/{sample}-tss-{norm}-{strand}.bw", norm=["spikenorm","libsizenorm"], sample=SAMPLES, strand=["SENSE","ANTISENSE","plus","minus"]),
        #datavis
        expand("datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bygroup.svg", annotation = config["annotations"], norm = ["spikenorm", "libsizenorm"], strand = ["SENSE", "ANTISENSE"]),
        #quality control
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all", "passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-tss-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status = ["all", "passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-tss-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status = ["all", "passing"]),
        #call DE bases/clusters 
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-{{type}}-qcplots-libsizenorm.svg", zip, condition=conditiongroups, control=controlgroups),type=["base", "cluster"]),
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-{{type}}-qcplots-spikenorm.svg", zip, condition=conditiongroups_si, control=controlgroups_si),type=["base", "cluster"]),
        #base and cluster distances
        # expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-distances-libsizenorm-{{direction}}.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"]),
        #call DE genic
        expand("diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-qcplots-libsizenorm.svg", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-qcplots-spikenorm.svg", zip, condition=conditiongroups_si, control=controlgroups_si),
        #
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-{{type}}-results-libsizenorm-{{direction}}.bed", zip, condition=conditiongroups, control=controlgroups),type=["base","cluster"], direction=["up","down"]),
        expand(expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-{{type}}-results-spikenorm-{{direction}}.bed", zip, condition=conditiongroups_si, control=controlgroups_si),type=["base","cluster"], direction=["up","down"]),
        expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-cluster-results-libsizenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups, control=controlgroups), direction = ["up","down"], category=CATEGORIES),
        expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-cluster-results-spikenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up","down"], category=CATEGORIES),
        #find intragenic ORFs
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intragenic-orfs/{condition}-v-{control}-libsizenorm-{{direction}}-intragenic-orfs.tsv", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"]),
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intragenic-orfs/{condition}-v-{control}-spikenorm-{{direction}}-intragenic-orfs.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"]),
        #MEME-ChIP
        # expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-spikenorm-{{direction}}-{{category}}-motifs/index.html", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"], category = CATEGORIES),
        # expand(expand("diff_exp/{condition}-v-{control}/{{category}}/{condition}-v-{control}-libsizenorm-{{direction}}-{{category}}-motifs/index.html", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"], category = CATEGORIES),
        # expand(expand("diff_exp/{condition}-v-{control}/{{type}}/{{type}}-v-genic/{condition}-v-{control}-{{type}}-v-genic-spikenorm.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), type=["antisense", "convergent", "divergent", "intragenic"]),
        # expand(expand("diff_exp/{condition}-v-{control}/{{type}}/{{type}}-v-genic/{condition}-v-{control}-{{type}}-v-genic-libsizenorm.tsv", zip, condition=conditiongroups, control=controlgroups), type=["antisense", "convergent", "divergent", "intragenic"]),
        # intrafreq
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-libsizenorm-{{direction}}-freqperORF.svg", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"]),
        # expand(expand("diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-spikenorm-{{direction}}-freqperORF.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"]),

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
        SIplmin = "coverage/counts/spikein/{sample}-tss-SI-counts-plmin.bedgraph"
    output:
        spikePlus = "coverage/spikenorm/{sample}-tss-spikenorm-plus.bedgraph",
        spikeMinus = "coverage/spikenorm/{sample}-tss-spikenorm-minus.bedgraph",
        libnormPlus = "coverage/libsizenorm/{sample}-tss-libsizenorm-plus.bedgraph",
        libnormMinus = "coverage/libsizenorm/{sample}-tss-libsizenorm-minus.bedgraph"
    params: scalefactor = config["spikein-pct"]
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
        SIplmin = "coverage/counts/spikein/{sample}-tss-SI-counts-plmin.bedgraph"
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
        samplelist = lambda wildcards : list({k:v for (k,v) in SAMPLES.items() if v["spikein"]=="y"}.keys()) if wildcards.status=="all" else list({k:v for (k,v) in PASSING.items() if v["spikein"]=="y"}.keys()),
        conditions = config["comparisons"]["spikenorm"]["conditions"],
        controls = config["comparisons"]["spikenorm"]["controls"],
    script: "scripts/plotsipct.R"

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
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus} > {output.antisense}) &>> {log}
        """

rule make_stranded_sicounts_bedgraph:
    input:
        plus = "coverage/counts/spikein/{sample}-tss-SI-counts-plus.bedgraph",
        minus = "coverage/counts/spikein/{sample}-tss-SI-counts-minus.bedgraph"
    output:
        sense = "coverage/counts/spikein/{sample}-tss-SI-counts-SENSE.bedgraph"
    log: "logs/make_stranded_sicounts_bedgraph/make_stranded_sicounts_bedgraph-{sample}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
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

rule bg_to_bw:
    input:
        bedgraph = "coverage/{norm}/{sample}-tss-{norm}-{strand}.bedgraph",
        chrsizes = lambda wildcards: config["genome"]["chrsizes"] if (wildcards.strand=="plus" or wildcards.strand=="minus") else os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"
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
        heatmap_sample = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bysample.svg",
        heatmap_group = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bygroup.svg",
    params:
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
    log: "logs/union_bedgraph-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {name_string} | bash scripts/cleanUnionbedg.sh | pigz > {output.exp}) &> {log}
        """

rule union_bedgraph_si_counts:
    input:
        si = expand("coverage/counts/spikein/{sample}-tss-SI-counts-SENSE.bedgraph", sample=SAMPLES),
    output:
        si = "coverage/counts/spikein/union-bedgraph-allsamples-si-counts.tsv.gz",
    params:
    log: "logs/union_bedgraph_si_counts.log"
    shell: """
        (bedtools unionbedg -i {input.si} -header -names {name_string} | bash scripts/cleanUnionbedg.sh | pigz > {output.si}) &> {log}
        """

def plotcorrsamples(wildcards):
    dd = SAMPLES if wildcards.status=="all" else PASSING
    if wildcards.condition=="all":
        if wildcards.norm=="libsizenorm": #condition==all,norm==lib
            return list(dd.keys())
        else: #condition==all,norm==spike
            return list({k:v for (k,v) in dd.items() if v["spikein"]=="y"}.keys())
    elif wildcards.norm=="libsizenorm": #condition!=all;norm==lib
        return list({k:v for (k,v) in dd.items() if v["group"]==wildcards.control or v["group"]==wildcards.condition}.keys())
    else: #condition!=all;norm==spike
        return list({k:v for (k,v) in dd.items() if (v["group"]==wildcards.control or v["group"]==wildcards.condition) and v["spikein"]=="y"}.keys())

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

#NOTE: need to check whether the median of ratios normalization is okay (e.g. for spt6 samples)
rule call_de_bases:
    input:
        counts = "coverage/counts/union-bedgraph-allsamples-counts.tsv.gz",
        sicounts = "coverage/counts/spikein/union-bedgraph-allsamples-si-counts.tsv.gz"
    params:
        samples = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}.keys()),
        groups = lambda wildcards : [PASSING[x]["group"] for x in {k:v for (k,v) in PASSING.items() if (v["group"]==wildcards.control or v["group"]==wildcards.condition)}],
        alpha = config["deseq"]["fdr"],
    output:
        results = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-results-{norm}-all.tsv",
        #need to write out norm counts here or just in the total qc?
        normcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-counts-sfnorm-{norm}.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-counts-rlog-{norm}.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-qcplots-{norm}.svg"
    script:
        "scripts/call_de_bases2.R"

rule separate_de_results:
    input:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-results-{norm}-all.tsv",
    output:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-results-{norm}-up.tsv",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-results-{norm}-down.tsv",
    params:
        alpha = config["deseq"]["fdr"]
    log: "logs/separate_de_results/separate_de_results-{condition}-v-{control}-{type}-{norm}.log"
    shell: """
        (awk -v fdr={params.alpha} -v uout={output.up} -v dout={output.down} 'BEGIN{{FS=OFS="\t"}} $7<fdr{{if ($3<0) {{print > dout}} else {{print > uout}} }}' {input}) &> {log}
        """

rule de_bases_to_bed:
    input:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-results-{norm}-{direction}.tsv",
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-results-{norm}-{direction}.bed",
    log: "logs/de_results_to_bed/de_results_to_bed-{condition}-v-{control}-base-{norm}-{direction}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1, -log($7)/log(10)}}' {input} | awk -v dir={wildcards.direction} -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1"-"$2, $3, $4, dir"_"NR, $5, "+"}} $2=="minus"{{print $1"-"$2, $3, $4, dir"_"NR, $5, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule merge_de_bases_to_clusters:
    input:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-results-{norm}-up.bed",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-results-{norm}-down.bed"
    output:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-allclusters-{norm}-up.bed",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-allclusters-{norm}-down.bed",
        combined = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-allclusters-{norm}-combined.bed"
    params:
        mergedist = config["cluster-merge-distance"]
    log: "logs/merge_de_bases_to_clusters/merge_de_bases_to_clusters-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools merge -s -d {params.mergedist} -i {input.up} | LC_COLLATE=C sort -k1,1 -k2,2n | tee {output.up} | cat - <(bedtools merge -s -d {params.mergedist} -i {input.down} | LC_COLLATE=C sort -k1,1 -k2,2n | tee {output.down}) | LC_COLLATE=C sort -k1,1 -k2,2n > {output.combined}) &> {log}
        """

rule map_counts_to_clusters:
    input:
        bed = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-allclusters-{norm}-combined.bed",
        bg = "coverage/counts/{sample}-tss-counts-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/{sample}-allclustercounts-{norm}.tsv")
    log: "logs/map_counts_to_clusters/map_counts_to_clusters-{condition}-v-{control}-{sample}-{norm}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-"$2"-"$3, $5}}' &> {output}) &> {log}
        """

rule get_cluster_counts:
    input:
        lambda wildcards : ["diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/" + x + "-allclustercounts-" + wildcards.norm + ".tsv" for x in list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)})]
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-allclusters-{norm}-counts.tsv"
    params:
        n = lambda wildcards: 2*len({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}),
        names = lambda wildcards: "\t".join(list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}.keys()))
    log: "logs/get_cluster_counts/get_cluster_counts-{condition}-v-{control}-{norm}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule call_de_clusters:
    input:
        clustercounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-allclusters-{norm}-counts.tsv",
        libcounts = lambda wildcards : "coverage/counts/union-bedgraph-allsamples-counts.tsv.gz" if wildcards.norm=="libsizenorm" else "coverage/counts/spikein/union-bedgraph-allsamples-si-counts.tsv.gz",
    params:
        samples = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}.keys()),
        groups = lambda wildcards : [PASSING[x]["group"] for x in {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}],
        alpha = config["deseq"]["fdr"]
    output:
        results = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-all.tsv",
        #need to write out norm counts here or just in the total qc?
        normcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-counts-sfnorm-{norm}.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-counts-rlog-{norm}.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-qcplots-{norm}.svg"
    script:
        "scripts/call_de_clusters2.R"

#NOTE: column 5 for down tables vs column 5 for up tables is to the negative sign in the fold-change
rule de_clusters_to_bed:
    input:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-up.tsv",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-down.tsv",
    output:
        up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-up.bed",
        down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-down.bed",
    log: "logs/de_results_to_bed/de_results_to_bed-{condition}-v-{control}-base-{norm}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $3":"(-log($7)/log(10))}}' {input.up} | awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1, $3, $4, "up_"NR, $5, "+"}} $2=="minus"{{print $1, $3, $4, "up_"NR, $5, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.up}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $3":"(-log($7)/log(10))}}' {input.down} | awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1, $3, $4, "down_"NR, $6, "+"}} $2=="minus"{{print $1, $3, $4, "down_"NR, $6, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.down}) &>> {log}
        """

#NOTE: the clusters are 'all' clusters, not just those that turn out to be DE (though this is most of them)
# rule extract_base_cluster_dist:
#     input:
#         bases = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-results-{norm}-{direction}.bed",
#         clusters = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-allclusters-{norm}-{direction}.bed",
#     output:
#         bases = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-distances-{norm}-{direction}.tsv",
#         clusters = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-distances-{norm}-{direction}.tsv",
#     log: "logs/extract_base_cluster_dist/extract_base_cluster_dist-{condition}-v-{control}-{norm}-{direction}.log"
#     shell: """
#         (bedtools closest -s -d -io -t first -a {input.bases} -b {input.bases} | cut -f13 > {output.bases}) &> {log}
#         (bedtools closest -d -io -t first -a {input.clusters} -b {input.clusters} | cut -f9> {output.clusters}) &>> {log}
#         """
# rule vis_base_cluster_dist:
#     input:
#         upbases = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-distances-{norm}-up.tsv",
#         dnbases = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-base-distances-{norm}-down.tsv",
#         upclusters = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-distances-{norm}-up.tsv",
#         dnclusters = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-distances-{norm}-down.tsv",
#     output:

rule map_counts_to_genic:
    input:
        bed = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        bg = "coverage/counts/{sample}-tss-counts-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/all_genic/{sample}-allgeniccounts.tsv")
    log: "logs/map_counts_to_genic/map_counts_to_genic-{condition}-v-{control}-{sample}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input.bed} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum | cut -f4,7 > {output}) &> {log}
        """

rule get_genic_counts:
    input:
        lambda wildcards : ["diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/all_genic/" + x + "-allgeniccounts.tsv" for x in list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)})]
    output:
        "diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-counts.tsv"
    params:
        n = lambda wildcards: 2*len({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}),
        names = lambda wildcards: "\t".join(list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}.keys()))
    log: "logs/get_genic_counts/get_genic_counts-{condition}-v-{control}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule call_de_genic:
    input:
        clustercounts = "diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-counts.tsv",
        libcounts = lambda wildcards : "coverage/counts/union-bedgraph-allsamples-counts.tsv.gz" if wildcards.norm=="libsizenorm" else "coverage/counts/spikein/union-bedgraph-allsamples-si-counts.tsv.gz",
    params:
        samples = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}.keys()),
        groups = lambda wildcards : [PASSING[x]["group"] for x in {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}],
        alpha = config["deseq"]["fdr"]
    output:
        results = "diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-results-{norm}-all.tsv",
        #need to write out norm counts here or just in the total qc?
        normcounts = "diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-counts-sfnorm-{norm}.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-counts-rlog-{norm}.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/all_genic/{condition}-v-{control}-allgenic-qcplots-{norm}.svg"
    script:
        "scripts/call_de_clusters2.R"

#TODO: add a cat statement to add a header for all 'class' tsv (make sure to check class to bed afterwards)
rule get_putative_intragenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-{direction}.bed",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-cluster-results-{norm}-{direction}-intragenic.tsv"
    log: "logs/get_putative_intragenic/get_putative_intragenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}}' | sort -k10,10nr > {output}) &> {log}
        """

rule get_intragenic_frequency:
    input:
        orfs = config["genome"]["orf-annotation"],
        intrabed = "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-de-clusters-{norm}-{direction}-intragenic.bed"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-de-clusters-{norm}-{direction}-intrafreq.tsv"
    shell: """
        bedtools intersect -a {input.orfs} -b {input.intrabed} -c -s > {output}
        """

rule plot_intragenic_frequency:
    input:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-de-clusters-{norm}-{direction}-intrafreq.tsv"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-intragenic-{norm}-{direction}-freqperORF.svg"
    log: "logs/plot_intragenic_frequency/plot_intragenic_frequency-{condition}-v-{control}-{norm}-{direction}.log"
    script: "scripts/intrafreq.R"


#(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\torf_start\torf_end\torf_name\tpeak_lfc\tpeak_significance\tdist_peak_to_atg\n$
rule get_putative_antisense:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-{direction}.bed",
        transcripts = config["genome"]["transcripts"]
    output:
        "diff_exp/{condition}-v-{control}/antisense/{condition}-v-{control}-cluster-results-{norm}-{direction}-antisense.tsv"
    log : "logs/get_putative_antisense/get_putative_antisense-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcripts} -wo -S | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}}' | sort -k10,10nr > {output}) &> {log}
        """
#(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance\tdist_peak_to_senseTSS\n$

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

rule get_putative_genic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-{direction}.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/genic/{condition}-v-{control}-cluster-results-{norm}-{direction}-genic.tsv"
    log : "logs/get_putative_genic/get_putative_genic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6}}' | sort -k10,10nr > {output}) &> {log}
        """
#(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance\n$

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
        (bedtools slop -s -l {params.genic_up} -r 0 -i {input.transcripts} -g <(sort -k1,1 {input.chrsizes})| sort -k1,1 -k2,2n | bedtools complement -i stdin -g <(sort -k1,1 {input.chrsizes}) > {output}) &> {log}
        """
        #(sort -k1,1 {input.chrsizes} > .chrsizes.temp) &> {log}
        #(rm .chrsizes.temp) &>> {log}

rule get_putative_intergenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-{direction}.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/intergenic/{condition}-v-{control}-cluster-results-{norm}-{direction}-intergenic.tsv"
    log : "logs/get_putative_intergenic/get_putative_intergenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo | awk 'BEGIN{{FS="\t|:";OFS="\t"}}{{print $1, $7, $2, $3, $4, $9, $10, ".", $5, $6}}'| sort -k10,10nr > {output}) &> {log}
        """
#(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\tregion_start\tregion_end\tregion_name\tpeak_lfc\tpeak_significance\n

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

rule get_putative_convergent:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-{direction}.bed",
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/convergent/{condition}-v-{control}-cluster-results-{norm}-{direction}-convergent.tsv"
    log : "logs/get_putative_convergent/get_putative_convergent-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}}' | sort -k10,10nr > {output}) &> {log}
        """
#(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance\tdist_peak_to_senseTSS\n$

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

rule get_putative_divergent:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-cluster-results-{norm}-{direction}.bed",
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_exp/{condition}-v-{control}/divergent/{condition}-v-{control}-cluster-results-{norm}-{direction}-divergent.tsv"
    log : "logs/get_putative_divergent/get_putative_divergent-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -wo -s | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $7=="+"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}} $7=="-"{{print $1, $7, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}}' | sort -k10,10nr > {output}) &> {log}
        """
#(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance\tdist_peak_to_senseTSS\n$

rule get_category_bed:
    input:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-cluster-results-{norm}-{direction}-{category}.tsv"
    output:
        "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-cluster-results-{norm}-{direction}-{category}.bed"
    log: "logs/get_category_bed/get_category_bed-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $3, $4, $5, $10, $2}}' {input} | sort -k1,1 -k2,2n  > {output}) &> {log}
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
