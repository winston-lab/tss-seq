#!/usr/bin/env python
import os

configfile: "config.yaml"

SAMPLES = config["samples"]
name_string = " ".join(SAMPLES)
PASSING = {k:v for (k,v) in SAMPLES.items() if v["pass-qc"] == "pass"}
pass_string = " ".join(PASSING)

controlgroups = config["comparisons"]["libsizenorm"]["controls"]
conditiongroups = config["comparisons"]["libsizenorm"]["conditions"]
controlgroups_si = config["comparisons"]["spikenorm"]["controls"]
conditiongroups_si = config["comparisons"]["spikenorm"]["conditions"]

localrules: all,
            make_stranded_genome,
            make_stranded_bedgraph,
            make_stranded_sicounts_bedgraph,
            make_stranded_annotations,
            make_bigwig_for_deeptools,
            gzip_deeptools_table,
            melt_matrix,
            cat_matrices,
            union_bedgraph,
            union_bedgraph_cond_v_ctrl,
            cat_strands,
            separate_de_bases,
            de_bases_to_bed,
            merge_de_bases_to_clusters,
            cat_cluster_strands,
            #union_cluster_bedgraph,
            get_cluster_counts,
            extract_base_distances,
            separate_de_clusters,
	    de_clusters_to_bed,
            map_counts_to_clusters,
            get_putative_intragenic

rule all:
    input:
        expand("qual_ctrl/fastqc/raw/{sample}", sample=SAMPLES),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip", sample=SAMPLES),
        expand("datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bygroup.png", annotation = config["annotations"], norm = ["spikenorm", "libsizenorm"], strand = ["SENSE", "ANTISENSE"]),
        "qual_ctrl/all/all-pca-scree-libsizenorm.png",
        "qual_ctrl/passing/passing-pca-scree-libsizenorm.png",
        #expand("diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-spikenorm.tsv", zip, condition=conditiongroups_si, control=controlgroups_si),
        #expand("diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-libsizenorm.tsv", zip, condition=conditiongroups, control=controlgroups),
        #expand("diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-spikenorm-down.bed", zip, condition=conditiongroups_si, control=controlgroups_si),
        #expand("diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-libsizenorm-down.bed", zip, condition=conditiongroups, control=controlgroups),
        expand(expand("diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-de-clusters-spikenorm-{{direction}}-intragenic.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["up", "down"]),
        expand(expand("diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-de-clusters-libsizenorm-{{direction}}-intragenic.tsv", zip, condition=conditiongroups, control=controlgroups), direction = ["up", "down"]),
       
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
    log: "logs/bowtie2_build.log"
    shell: """
        (bowtie2-build {input.fasta} ../genome/bowtie2_indexes/{params.name}) &> {log}
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
        SIplmin = temp("coverage/counts/spikein/{sample}-tss-SI-counts-plmin.bedgraph"),
        SIpl = "coverage/counts/spikein/{sample}-tss-SI-counts-plus.bedgraph",
        SImin = "coverage/counts/spikein/{sample}-tss-SI-counts-minus.bedgraph",
        plmin = temp("coverage/counts/{sample}-tss-counts-plmin.bedgraph"),
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
        (rm coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-*.tmp) &>> {log} 
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2, $3, $4}}' {input.minus} > coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-plus.tmp) &>> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-minus", $2, $3, $4}}' {input.plus} > coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-minus.tmp) &>> {log}
        (cat coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-plus.tmp coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-minus.tmp | LC_COLLATE=C sort -k1,1 -k2,2n > {output.antisense}) &>> {log} 
        rm coverage/{wildcards.norm}/{wildcards.sample}-{wildcards.norm}-*.tmp 
        """

rule make_stranded_sicounts_bedgraph:
    input:
        plus = "coverage/counts/spikein/{sample}-tss-SI-counts-plus.bedgraph",
        minus = "coverage/counts/spikein/{sample}-tss-SI-counts-minus.bedgraph"
    output:
        sense = "coverage/counts/spikein/{sample}-tss-SI-counts-SENSE.bedgraph"
    log: "logs/make_stranded_sicounts_bedgraph/make_stranded_sicounts_bedgraph-{sample}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2, $3, $4}}' {input.plus} > coverage/counts/spikein/{wildcards.sample}-counts-plus.tmp) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-minus", $2, $3, $4}}' {input.minus} > coverage/counts/spikein/{wildcards.sample}-counts-minus.tmp) &>> {log}
        (cat coverage/counts/spikein/{wildcards.sample}-counts-plus.tmp coverage/counts/spikein/{wildcards.sample}-counts-minus.tmp | LC_COLLATE=C sort -k1,1 -k2,2n > {output.sense}) &>> {log} 
        (rm coverage/counts/spikein/{wildcards.sample}-counts-*.tmp) &>> {log}
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

rule gzip_deeptools_table:
    input:
        tsv = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv",
        mat = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.mat.gz"
    output:
        "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{strand}.tsv.gz"
    log: "logs/gzip_deeptools_table/gzip_deeptools_table-{annotation}-{sample}-{norm}-{strand}.log"
    shell: """
        (pigz -f {input.tsv}) &> {log}
        (rm {input.mat}) &>> {log}
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
        heatmap_sample = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bysample.png",
        heatmap_group = "datavis/{annotation}/{norm}/tss-{annotation}-{norm}-{strand}-heatmap-bygroup.png",
    params:
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
        #metagene_palette = lambda wildcards : config["annotations"][wildcards.annotation]["metagene_palette"],
        #avg_heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["avg_heatmap_cmap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plotHeatmaps.R"

rule union_bedgraph:
    input:
        exp = expand("coverage/counts/{sample}-tss-counts-SENSE.bedgraph", sample=SAMPLES),
        si = expand("coverage/counts/spikein/{sample}-tss-SI-counts-SENSE.bedgraph", sample=SAMPLES),
        pass_exp = expand("coverage/counts/{sample}-tss-counts-SENSE.bedgraph", sample=PASSING),
        pass_si = expand("coverage/counts/spikein/{sample}-tss-SI-counts-SENSE.bedgraph", sample=PASSING),
    output:
        exp = "coverage/counts/union-bedgraph-allsamples.txt",    
        si = "coverage/counts/spikein/union-bedgraph-si-allsamples.txt",
        pass_exp = "coverage/counts/union-bedgraph-passing.txt",    
        pass_si = "coverage/counts/spikein/union-bedgraph-si-passing.txt"
    params:
        allminreads = config["minreads"]*len(SAMPLES),
        si_allminreads = config["minreads"]*len(SAMPLES)/10,
        passminreads = config["minreads"]*len(PASSING),
        si_passminreads = config["minreads"]*len(PASSING)/10
    log: "logs/union_bedgraph.log"
    #TODO: write this into a script
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {name_string} |
        awk -v awkmin={params.allminreads} 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0}} {{t=0; for(i=4; i<=NF; i++) t+=$i}} t>awkmin{{print $0}}' > .union-bedgraph-allsamples.temp) &> {log}
        (cut -f1-3 .union-bedgraph-allsamples.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}'  > .positions.txt) &>> {log}
        (cut -f4- .union-bedgraph-allsamples.temp > .values.txt) &>> {log}
        (paste .positions.txt .values.txt > {output.exp}) &>> {log}
        (rm .union-bedgraph-allsamples.temp .positions.txt .values.txt) &>> {log}

        (bedtools unionbedg -i {input.si} -header -names {name_string} |
        awk -v awkmin={params.si_allminreads} 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0}} {{t=0; for(i=4; i<=NF; i++) t+=$i}} t>awkmin{{print $0}}' > .union-bedgraph-si-allsamples.temp) &>> {log}
        (cut -f1-3 .union-bedgraph-si-allsamples.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}'  > .positions.txt) &>> {log}
        (cut -f4- .union-bedgraph-si-allsamples.temp > .values.txt) &>> {log}
        (paste .positions.txt .values.txt > {output.si}) &>> {log}
        (rm .union-bedgraph-si-allsamples.temp .positions.txt .values.txt) &>> {log}

        (bedtools unionbedg -i {input.pass_exp} -header -names {pass_string} |
        awk -v awkmin={params.passminreads} 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0}} {{t=0; for(i=4; i<=NF; i++) t+=$i}} t>awkmin{{print $0}}' > .union-bedgraph-passing.temp) &>> {log}
        (cut -f1-3 .union-bedgraph-passing.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}'  > .positions.txt) &>> {log}
        (cut -f4- .union-bedgraph-passing.temp > .values.txt) &>> {log}
        (paste .positions.txt .values.txt > {output.pass_exp}) &>> {log}
        (rm .union-bedgraph-passing.temp .positions.txt .values.txt) &>> {log}

        (bedtools unionbedg -i {input.pass_si} -header -names {pass_string} |
        awk -v awkmin={params.si_passminreads} 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0}} {{t=0; for(i=4; i<=NF; i++) t+=$i}} t>awkmin{{print $0}}' > .union-bedgraph-si-passing.temp) &>> {log}
        (cut -f1-3 .union-bedgraph-si-passing.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}'  > .positions.txt) &>> {log}
        (cut -f4- .union-bedgraph-si-passing.temp > .values.txt) &>> {log}
        (paste .positions.txt .values.txt > {output.pass_si}) &>> {log}
        (rm .union-bedgraph-si-passing.temp .positions.txt .values.txt) &>> {log}
        """

rule union_bedgraph_cond_v_ctrl:
    input:
        exp = lambda wildcards : expand("coverage/counts/{sample}-tss-counts-SENSE.bedgraph", sample = {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}),
        si = lambda wildcards : expand("coverage/counts/spikein/{sample}-tss-SI-counts-SENSE.bedgraph", sample = {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)})
    output:
        exp = "coverage/counts/union-bedgraph-{condition}-v-{control}.txt",
        si = "coverage/counts/spikein/union-bedgraph-si-{condition}-v-{control}.txt"
    params:
        names = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}.keys()),
        minreads = lambda wildcards : config["minreads"]*len({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}),
        si_minreads = lambda wildcards : config["minreads"]*len({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)})/10
    log: "logs/union_bedgraph_cond_v_ctrl/union_bedgraph_{condition}-v-{control}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} |
        awk -v awkmin={params.minreads} 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0}} {{t=0; for(i=4; i<=NF; i++) t+=$i}} t>awkmin{{print $0}}' > .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp) &> {log}
        (cut -f1-3 .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}'  > .{wildcards.condition}-v-{wildcards.control}-positions.txt) &>> {log}
        (cut -f4- .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp > .{wildcards.condition}-v-{wildcards.control}-values.txt) &>> {log}
        (paste .{wildcards.condition}-v-{wildcards.control}-positions.txt .{wildcards.condition}-v-{wildcards.control}-values.txt > {output.exp}) &>> {log}
        (rm .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp .{wildcards.condition}-v-{wildcards.control}-positions.txt .{wildcards.condition}-v-{wildcards.control}-values.txt) &>> {log}

        (bedtools unionbedg -i {input.si} -header -names {params.names} |
        awk -v awkmin={params.si_minreads} 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0}} {{t=0; for(i=4; i<=NF; i++) t+=$i}} t>awkmin{{print $0}}' > .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp) &>> {log}
        (cut -f1-3 .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}'  > .{wildcards.condition}-v-{wildcards.control}-positions.txt) &>> {log}
        (cut -f4- .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp > .{wildcards.condition}-v-{wildcards.control}-values.txt) &>> {log}
        (paste .{wildcards.condition}-v-{wildcards.control}-positions.txt .{wildcards.condition}-v-{wildcards.control}-values.txt > {output.si}) &>> {log}
        (rm .unionbedg-{wildcards.condition}-v-{wildcards.control}.temp .{wildcards.condition}-v-{wildcards.control}-positions.txt .{wildcards.condition}-v-{wildcards.control}-values.txt) &>> {log}
        """

rule deseq_initial_qc:
    input:
        exp = "coverage/counts/union-bedgraph-allsamples.txt",
        si = "coverage/counts/spikein/union-bedgraph-si-allsamples.txt"
    params:
        alpha = config["deseq"]["fdr"],
        samples = list(SAMPLES.keys()),
        samplegroups = [SAMPLES[x]["group"] for x in SAMPLES],
        nospikein = list({k:v for (k,v) in SAMPLES.items() if (v["spikein"] == "n")}.keys())
    output:
        exp_size_v_sf = protected("qual_ctrl/all/all-libsize-v-sizefactor-experimental.png"),
        si_size_v_sf = protected("qual_ctrl/all/all-libsize-v-sizefactor-spikein.png"),
        si_pct = protected("qual_ctrl/all/all-spikein-pct.png"),
        corrplot_spikenorm = protected("qual_ctrl/all/all-pairwise-correlation-spikenorm.png"),
        corrplot_libsizenorm = protected("qual_ctrl/all/all-pairwise-correlation-libsizenorm.png"),
        #count_heatmap_spikenorm = protected("qual_ctrl/all/all-de-bases-heatmap-spikenorm.png"),
        #count_heatmap_libsizenorm = protected("qual_ctrl/all/all-de-bases-heatmap-libsizenorm.png"),
        dist_heatmap_spikenorm = protected("qual_ctrl/all/all-sample-dists-spikenorm.png"),
        dist_heatmap_libsizenorm = protected("qual_ctrl/all/all-sample-dists-libsizenorm.png"),
        pca_spikenorm = protected("qual_ctrl/all/all-pca-spikenorm.png"),
        scree_spikenorm = protected("qual_ctrl/all/all-pca-scree-spikenorm.png"),
        pca_libsizenorm = protected("qual_ctrl/all/all-pca-libsizenorm.png"),
        scree_libsizenorm = protected("qual_ctrl/all/all-pca-scree-libsizenorm.png")
    script:
       "scripts/deseq2_qc.R"

rule deseq_passing_qc:
    input:
        exp = "coverage/counts/union-bedgraph-passing.txt",
        si = "coverage/counts/spikein/union-bedgraph-si-passing.txt"
    params:
        alpha = config["deseq"]["fdr"],
        samples = list(PASSING.keys()),
        samplegroups = [PASSING[x]["group"] for x in PASSING],
        nospikein = list({k:v for (k,v) in PASSING.items() if v["spikein"] == "n"}.keys())
    output:
        exp_size_v_sf = protected("qual_ctrl/passing/passing-libsize-v-sizefactor-experimental.png"),
        si_size_v_sf = protected("qual_ctrl/passing/passing-libsize-v-sizefactor-spikein.png"),
        si_pct = protected("qual_ctrl/passing/passing-spikein-pct.png"),
        corrplot_spikenorm = protected("qual_ctrl/passing/passing-pairwise-correlation-spikenorm.png"),
        corrplot_libsizenorm = protected("qual_ctrl/passing/passing-pairwise-correlation-libsizenorm.png"),
        #count_heatmap_spikenorm = protected("qual_ctrl/passing/passing-de-bases-heatmap-spikenorm.png"),
        #count_heatmap_libsizenorm = protected("qual_ctrl/passing/passing-de-bases-heatmap-libsizenorm.png"),
        dist_heatmap_spikenorm = protected("qual_ctrl/passing/passing-sample-dists-spikenorm.png"),
        dist_heatmap_libsizenorm = protected("qual_ctrl/passing/passing-sample-dists-libsizenorm.png"),
        pca_spikenorm = protected("qual_ctrl/passing/passing-pca-spikenorm.png"),
        scree_spikenorm = protected("qual_ctrl/passing/passing-pca-scree-spikenorm.png"),
        pca_libsizenorm = protected("qual_ctrl/passing/passing-pca-libsizenorm.png"),
        scree_libsizenorm = protected("qual_ctrl/passing/passing-pca-scree-libsizenorm.png")
    script:
        "scripts/deseq2_qc.R"

rule call_de_bases_cond_v_ctrl:
    input:
        exp = "coverage/counts/union-bedgraph-{condition}-v-{control}.txt",
        si = "coverage/counts/spikein/union-bedgraph-si-{condition}-v-{control}.txt"
    params:
        alpha = config["deseq"]["fdr"],
        samples = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}.keys()),
        samplegroups = lambda wildcards : [PASSING[x]["group"] for x in {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}],
        nospikein = lambda wildcards : list({k:v for (k,v) in PASSING.items() if ((v["group"]== wildcards.control or v["group"]==wildcards.condition) and v["spikein"] == "n")}.keys())
    output:
        #exp_size_v_sf = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-libsize-v-sizefactor-experimental.png",
        #si_size_v_sf = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-libsize-v-sizefactor-spikein.png",
        #si_pct = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-spikein-pct.png",
        corrplot_spikenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-pairwise-correlation-spikenorm.png",
        corrplot_libsizenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-pairwise-correlation-libsizenorm.png",
        #count_heatmap_spikenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-de-bases-heatmap-spikenorm.png",
        #count_heatmap_libsizenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-de-bases-heatmap-libsizenorm.png",
        dist_heatmap_spikenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-sample-dists-spikenorm.png",
        dist_heatmap_libsizenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-sample-dists-libsizenorm.png",
        pca_spikenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-pca-spikenorm.png",
        scree_spikenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-pca-scree-spikenorm.png",
        pca_libsizenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-pca-libsizenorm.png",
        scree_libsizenorm = "qual_ctrl/{condition}-v-{control}/{condition}-v-{control}-pca-scree-libsizenorm.png",
        de_spikenorm_path = protected("diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-spikenorm.tsv"),
        de_libsizenorm_path = protected("diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-libsizenorm.tsv")
    script:
        "scripts/call_de_bases.R"

rule separate_de_bases:
    input:
        "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}.tsv"
    output:
        up = "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}-up.tsv",
        down = "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}-down.tsv"
    log: "logs/separate_de_bases/separate_de_bases-{condition}-v-{control}-{norm}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} $3>=0 {{print $0}}' {input} > {output.up}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}} $3<0 {{print $0}}' {input} > {output.down}) &>> {log}
        """

rule de_bases_to_bed:
    input:
        up = "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}-up.tsv",
        down = "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}-down.tsv"
    output:
        up = "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}-up.bed", 
        down = "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}-down.bed" 
    log: "logs/de_bases_to_bed/de_bases_to_bed-{condition}-v-{control}-{norm}.log"
    shell: """
        (tail -n +2 {input.up} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, -log($7)/log(10)}}' | awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1"-"$2, $3, $4, "up_"NR, $5, "+"}} $2=="minus"{{print $1"-"$2, $3, $4, "up_"NR, $5, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.up}) &> {log}
        (tail -n +2 {input.down} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, -log($7)/log(10)}}'| awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1"-"$2, $3, $4, "down_"NR, $5, "+"}} $2=="minus"{{print $1"-"$2, $3, $4, "down_"NR, $5, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.down}) &>> {log}
        """

#TODO: rewrite this to extract base and cluster distances
#rule extract_base_distances:
#    input:
#        up = "diff_exp/de_bases/de-bases-{norm}-up.bed",
#        down = "diff_exp/de_bases/de-bases-{norm}-down.bed"
#    output:
#        "diff_exp/de_bases/base-distances-{norm}.tsv"
#    shell: """
#        bedtools closest -s -d -io -t first -a {input.up} -b {input.up} | cut -f13 > diff_exp/.{wildcards.norm}-basedistances-up.temp
#        bedtools closest -s -d -io -t first -a {input.down} -b {input.down} | cut -f13 > diff_exp/.{wildcards.norm}-basedistances-down.temp
#        cat diff_exp/.{wildcards.norm}-basedistances-up.temp diff_exp/.{wildcards.norm}-basedistances-down.temp > {output}
#        rm diff_exp/.{wildcards.norm}*.temp
#        """

#rule vis_base_distances

rule merge_de_bases_to_clusters:
    input:
        "diff_exp/{condition}-v-{control}/de_bases/{condition}-v-{control}-de-bases-{norm}-{direction}.bed"
    output:
        "diff_exp/{condition}-v-{control}/de_bases/allclusters/{condition}-v-{control}-allclusters-{norm}-{direction}.bed"
    params:
        mergedist = config["cluster-merge-distance"]
    log: "logs/merge_de_bases_to_clusters/merge_de_bases_to_clusters-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools merge -s -d {params.mergedist} -i {input} | LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule cat_cluster_strands:
    input:
        expand("diff_exp/{{condition}}-v-{{control}}/de_bases/allclusters/{{condition}}-v-{{control}}-allclusters-{{norm}}-{direction}.bed", direction = ["up", "down"]) 
    output:
        "diff_exp/{condition}-v-{control}/de_bases/allclusters/{condition}-v-{control}-allclusters-{norm}-combined.bed"
    log: "logs/cat_cluster_strands/cat_cluster_strands-{condition}-v-{control}-{norm}.log"
    shell: """
        (cat {input} | LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule map_counts_to_clusters:
    input:
        bed = "diff_exp/{condition}-v-{control}/de_bases/allclusters/{condition}-v-{control}-allclusters-{norm}-combined.bed",
        bg = "coverage/counts/{sample}-tss-counts-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/de_clusters/{sample}-allclusters-{norm}.tsv")
    log: "logs/map_counts_to_clusters/map_counts_to_clusters-{condition}-v-{control}-{sample}-{norm}.log"
    #shell: """
    #    (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | cut -f1,2,3,5 > {output}) &> {log}
    #    """
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-"$2"-"$3, $5}}' &> {output}) &> {log}
        """

rule get_cluster_counts:
    input:
        lambda wildcards : ["diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/de_clusters/" + x + "-allclusters-" + wildcards.norm + ".tsv" for x in list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)})]
    output:
        "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-{norm}-cluster-counts.tsv"
    shell: """
        bash scripts/recursivejoin.sh {input} > {output}
        """

#rule union_cluster_bedgraph:
#    input:
#        lambda wildcards : ["diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/de_clusters/" + x + "-allclusters-" + wildcards.norm + ".bedgraph" for x in list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)})]
#    output:
#        "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-union-bedgraph-clusters-{norm}.txt"
#    log: "logs/union_cluster_bedgraph/union_cluster_bedgraph-{condition}-v-{control}-{norm}.log"
#    shell: """
#        (bedtools unionbedg -i {input} > diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}.temp) &> {log}
#        (awk 'BEGIN{{FS="-|\t"; OFS=":"}} {{print $2, $1, $3, $4}}' diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}.temp > diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-base.temp) &>> {log}
#        (cut -f4- diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}.temp > diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-values.temp) &>> {log}
#        (paste diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-base.temp diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}-values.temp > {output}) &>> {log}
#        (rm diff_exp/{wildcards.condition}-v-{wildcards.control}/de_clusters/{wildcards.condition}-v-{wildcards.control}-{wildcards.norm}*temp) &>> {log}
#        """

rule call_de_clusters_spikenorm:
   input:
        clustercounts= "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-spikenorm-cluster-counts.tsv",
        libcounts = "coverage/counts/spikein/union-bedgraph-si-{condition}-v-{control}.txt"
   params:
        alpha = config["deseq"]["fdr"],
        samples = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}.keys()),
        samplegroups = lambda wildcards : [PASSING[x]["group"] for x in {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}]
   output:
        corrplot= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-pairwise-correlation-spikenorm.png",
        count_heatmap= "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-heatmap-spikenorm.png",
        dist_heatmap= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-sample-dists-spikenorm.png",
        pca= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-pca-spikenorm.png",
        scree= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-pca-scree-spikenorm.png",
        de_path = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-spikenorm.tsv",
   script:
        "scripts/call_de_clusters.R"

rule call_de_clusters_libsizenorm:
   input:
        clustercounts= "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-libsizenorm-cluster-counts.tsv",
        libcounts = "coverage/counts/union-bedgraph-{condition}-v-{control}.txt"
   params:
        alpha = config["deseq"]["fdr"],
        samples = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}.keys()),
        samplegroups = lambda wildcards : [PASSING[x]["group"] for x in {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]== wildcards.condition)}]
   output:
        corrplot= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-pairwise-correlation-libsizenorm.png",
        count_heatmap= "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-heatmap-libsizenorm.png",
        dist_heatmap= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-sample-dists-libsizenorm.png",
        pca= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-pca-libsizenorm.png",
        scree= "qual_ctrl/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-pca-scree-libsizenorm.png",
        de_path = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-libsizenorm.tsv",
   script:
        "scripts/call_de_clusters.R"

rule separate_de_clusters:
    input:
        "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}.tsv"
    output:
        up = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}-up.tsv",
        down = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}-down.tsv"
    log: "logs/separate_de_clusters/separate_de_clusters-{condition}-v-{control}-{norm}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} $3>=0 {{print $0}}' {input} > {output.up}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}} $3<0 {{print $0}}' {input} > {output.down}) &>> {log}
        """

rule de_clusters_to_bed:
    input:
        up = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}-up.tsv",
        down = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}-down.tsv"
    output:
       up = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}-up.bed", 
       down= "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}-down.bed" 
    log: "logs/de_clusters_to_bed/de_clusters_to_bed-{condition}-v-{control}-{norm}.log"
    shell: """
        (tail -n +2 {input.up} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, -log($7)/log(10)}}' | awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1, $3, $4, "up_"NR, $5, "+"}} $2=="minus"{{print $1, $3, $4, "up_"NR, $5, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.up}) &> {log}
        (tail -n +2 {input.down} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, -log($7)/log(10)}}'| awk -F '[-\t]' 'BEGIN{{OFS="\t"}} $2=="plus"{{print $1, $3, $4, "down_"NR, $5, "+"}} $2=="minus"{{print $1, $3, $4, "down_"NR, $5, "-"}}' | LC_COLLATE=C sort -k1,1 -k2,2n > {output.down}) &>> {log}
        """

rule get_putative_intragenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/de_clusters/{condition}-v-{control}-de-clusters-{norm}-{direction}.bed",
        orfs = config["genome"]["orf-annotation"] 
    output:
        "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-de-clusters-{norm}-{direction}-intragenic.tsv"
    log: "logs/get_putative_intragenic/get_putative_intragenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.orfs} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $6, $2, $3, $4, $8, $9, $10, $5, ((($2+1)+$3)/2)-$8}} $6=="-"{{print $1, $6, $2, $3, $4, $8, $9, $10, $5, $9-((($2+1)+$3)/2)}}' | sort -k9,9nr > {output}) &> {log}
        """
