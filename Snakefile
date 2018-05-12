#!/usr/bin/env python
import os
import subprocess
import itertools
from math import log2, log10

configfile: "config.yaml"

SAMPLES = config["samples"]
sisamples = {k:v for k,v in SAMPLES.items() if v["spikein"]=="y"}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]=="pass"}
sipassing = {k:v for k,v in PASSING.items() if v["spikein"]=="y"}

#groups which have at least two passing samples, so that they are valid for peakcalling and diff exp
validgroups = set([z for z in [PASSING[x]['group'] for x in PASSING] if [PASSING[x]['group'] for x in PASSING].count(z)>=2])
validgroups_si = set([z for z in [PASSING[x]['group'] for x in PASSING if PASSING[x]['spikein']=="y"] if [PASSING[x]['group'] for x in PASSING].count(z)>=2])
controlgroups = [v for k,v in config["comparisons"]["libsizenorm"].items() if v in validgroups]
conditiongroups = [k for k,v in config["comparisons"]["libsizenorm"].items() if k in validgroups]
controlgroups_si = [v for k,v in config["comparisons"]["spikenorm"].items() if v in validgroups_si]
conditiongroups_si = [k for k,v in config["comparisons"]["spikenorm"].items() if k in validgroups_si]

CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

FIGURES = config["figures"]

localrules:
    all,
    fastqc_aggregate,
    bowtie2_build,
    get_si_pct, plot_si_pct,
    make_stranded_genome, make_stranded_annotations,
    build_genic_annotation, build_convergent_annotation, build_divergent_annotation, build_intergenic_annotation,
    classify_peaks_genic, classify_peaks_intragenic, classify_peaks_antisense,
    classify_peaks_convergent, classify_peaks_divergent, classify_peaks_intergenic,
    combine_tss_peaks, map_counts_to_peaks, get_peak_counts,
    separate_de_peaks, de_peaks_to_bed,
    get_de_genic, get_de_intragenic, get_de_antisense,
    get_de_convergent, get_de_divergent, get_de_intergenic,
    cat_matrices,
    # get_de_intragenic_frequency
    # plot_de_intragenic_frequency
    # get_intra_orfs
    # separate_sig_de, get_de_category_bed,
    get_meme_de_peak_sequences,
    plot_seqlogos,
    # class_v_genic

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule all:
    input:
        #FastQC
        # 'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=SAMPLES),
        ##quality controls
        #"qual_ctrl/read_processing-loss.svg",
        #expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all", "passing"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-tss-{{status}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status = ["all", "passing"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-tss-{{status}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status = ["all", "passing"]),
        ##coverage
        #expand("coverage/{norm}/{sample}-tss-{norm}-{strand}.bw", norm=["spikenorm","libsizenorm", "counts", "sicounts"], sample=SAMPLES, strand=["SENSE","ANTISENSE","plus","minus"]),
        ##peakcalling on all samples
        #expand("peakcalling/sample_peaks/{sample}-exp-allpeaks.narrowPeak", sample=SAMPLES),
        #expand("peakcalling/sample_peaks/{sample}-si-allpeaks.narrowPeak", sample=sisamples),
        ##IDR for all groups which have at least two passing samples
        #expand("peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak", group=validgroups),
        #expand("peakcalling/{group}/{group}-si-idrpeaks.narrowPeak", group=validgroups_si),
        ##classify peaks into categories
        #expand("peakcalling/{group}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
        #expand("peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peakdistances.svg", zip, condition=conditiongroups + ["all"], control=controlgroups + ["all"]),
        ##combine called peaks for conditions vs control
        #expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-peaks.bed", zip, condition=conditiongroups, control=controlgroups),
        #expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}-si-peaks.bed", zip, condition=conditiongroups_si, control=controlgroups_si),
        ##differential expression of peaks
        #expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-experimental-peak-counts.tsv", zip, condition=conditiongroups, control=controlgroups),
        #expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-spikein-peak-counts.tsv", zip, condition=conditiongroups_si, control=controlgroups_si),
        #expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-results-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups),direction=["all", "up", "unchanged", "down"]),
        #expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{condition}-v-{control}_tss-seq-spikenorm-diffexp-results-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si),direction=["all", "up", "unchanged", "down"]),
        ###categorize differentially expressed peaks
        #expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-results-{{category}}-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups), direction = ["all","up","unchanged","down"], category="genic"),
        ## expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}-results-spikenorm-{{direction}}-{{category}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["all","up","unchanged","down"], category=CATEGORIES),
        ##differential expression summary
        #expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{condition}-v-{control}-libsizenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups, control=controlgroups), plot = ["summary", "maplot", "volcano"]),
        #expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{condition}-v-{control}-spikenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), plot = ["summary", "maplot", "volcano"]),
        ##distances of differentially expressed peaks
        #expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{{ttype}}/{condition}-v-{control}-relative-distances-libsizenorm-{{direction}}-{{ttype}}.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], ttype=["intragenic", "antisense"]),
        #expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{{ttype}}/{condition}-v-{control}-relative-distances-spikenorm-{{direction}}-{{ttype}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], ttype=["intragenic", "antisense"]),
        ##datavis
        #expand("datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{norm}-{strand}.tsv.gz", figure=FIGURES, norm=["spikenorm", "libsizenorm"], strand=["SENSE","ANTISENSE"]) if config["plot_figures"] else [],
        #expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/tss-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup-sense.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, status=["all","passing"]) if config["plot_figures"] else [],
        #expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/tss-{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup-sense.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, status=["all","passing"]) if config["plot_figures"] else [],
        ##correlations of transcript classes with genic TSSs
        #expand("diff_exp/{condition}-v-{control}/libsizenorm/genic_v_class/{condition}-v-{control}-libsizenorm-genic-v-class.svg", zip, condition=conditiongroups, control=controlgroups),
        #expand("diff_exp/{condition}-v-{control}/spikenorm/genic_v_class/{condition}-v-{control}-spikenorm-genic-v-class.svg", zip, condition=conditiongroups_si, control=controlgroups_si),
        ##enrichment of known motifs
        #"motifs/allmotifs.bed" if config["motifs"]["run_motif_analyses"] else [],
        #expand(expand("motifs/{condition}-v-{control}/libsizenorm/{{negative}}/{condition}-v-{control}_libsizenorm-{{direction}}-v-{{negative}}-{{category}}-motif_enrichment.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES) if config["motifs"]["run_motif_analyses"] else [],
        #expand(expand("motifs/{condition}-v-{control}/spikenorm/{{negative}}/{condition}-v-{control}_spikenorm-{{direction}}-v-{{negative}}-{{category}}-motif_enrichment.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES) if config["motifs"]["run_motif_analyses"] else [],
        ## MEME-ChIP
        #expand(expand("motifs/meme/{condition}-v-{control}/libsizenorm/{{region}}/{condition}-v-{control}-results-libsizenorm-{{direction}}-{{category}}-{{region}}-meme_chip/meme-chip.html", zip, condition=conditiongroups, control=controlgroups), region=["upstream","peak"], direction=["up", "down"], category=CATEGORIES) if config["motifs"]["meme-chip"]["run-meme-chip"] else [],
        #expand(expand("motifs/meme/{condition}-v-{control}/spikenorm/{{region}}/{condition}-v-{control}-results-spikenorm-{{direction}}-{{category}}-{{region}}-meme_chip/meme-chip.html", zip, condition=conditiongroups_si, control=controlgroups_si), region=["upstream","peak"], direction=["up", "down"], category=CATEGORIES) if config["motifs"]["meme-chip"]["run-meme-chip"] else [],
        ##GC pct coverage file
        #os.path.splitext(config["genome"]["fasta"])[0] + "-GC_pct.bw",
        ##gene ontology
        #expand(expand("gene_ontology/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}-GO-libsizenorm-{{direction}}-{{category}}-enriched-all.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down", "unchanged"], category=["genic", "intragenic", "antisense", "convergent", "divergent"]) if config["run_gene_ontology"] else [],
        #expand(expand("gene_ontology/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}-GO-spikenorm-{{direction}}-{{category}}-enriched-all.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down", "unchanged"], category=["genic", "intragenic", "antisense", "convergent", "divergent"]) if config["run_gene_ontology"] else [],
        ##sequence logos
        #expand("seq_logos/{group}/{group}-seqlogos.svg", group=set([PASSING[x]['group'] for x in PASSING]))
        ##find intragenic ORFs
        ##intragenic frequency per ORF

def plotcorrsamples(wc):
    if wc.condition=="all":
        if wc.norm=="libsizenorm": #condition==all,norm==lib
            return list(SAMPLES.keys())
        else: #condition==all,norm==spike
            return list(sisamples.keys())
    elif wc.norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in PASSING.items() if v["group"] in [wc.condition, wc.control]]
    else: #condition!=all;norm==spike
        return [k for k,v in sipassing.items() if v["group"] in [wc.control, wc.condition]]

def cluster_samples(status, norm, cluster_groups, cluster_strands):
    ll = []
    dd = SAMPLES if status=="all" else PASSING
    for group, strand in zip(cluster_groups, cluster_strands):
        sublist = [k for k,v in dd.items() if v["group"]==group] if norm=="libsizenorm" else [k for k,v in dd.items() if v["group"]==group and v["spikein"]=="y"]
        if strand in ["sense", "both"]:
            ll.append([sample + "-" + "sense" for sample in sublist])
        if strand in ["antisense", "both"]:
            ll.append([sample + "-" + "antisense" for sample in sublist])
    return(list(itertools.chain(*ll)))

rule fastqc_raw:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/raw/{sample}/{fname}/fastqc_data.txt"
    threads: config["threads"]
    log: "logs/fastqc/raw/fastqc_raw-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/raw/{wildcards.sample}) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw/{wildcards.sample} {input}) &>> {log}
        """

#in this order: remove adapter, remove 3' molecular barcode, do NextSeq quality trimming
#reads shorter than 6nt after cleaning are discarded
rule clean_reads:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    output:
        # fq = temp("fastq/cleaned/{sample}_tss-seq-trimmed.fastq.gz"),
        fq = "fastq/cleaned/{sample}_tss-seq-trimmed.fastq.gz",
        adapter = "logs/clean_reads/remove_adapter-{sample}.log",
        qual_trim = "logs/clean_reads/remove_3p_bc_and_trim-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"]
    threads: config["threads"]
    shell: """
        cutadapt --adapter={params.adapter} --minimum-length=18 {input} --cores={threads} 2> {output.adapter} | cutadapt --cut=-6 --nextseq-trim={params.trim_qual} --minimum-length=12 --output={output.fq} --cores={threads} - &> {output.qual_trim}
        """

rule remove_molecular_barcode:
    input:
        "fastq/cleaned/{sample}_tss-seq-trimmed.fastq.gz"
    output:
        fq = "fastq/cleaned/{sample}_tss-seq-clean.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    threads: config["threads"]
    log: "logs/remove_molecular_barcode/remove_molecular_barcode-{sample}.log"
    shell: """
        (python scripts/extractMolecularBarcode.py {input} {output.fq} {output.barcodes} {output.ligation}) &> {log}
        """
        # (python scripts/extractMolecularBarcode.py {input} fastq/cleaned/{wildcards.sample}_tss-seq-clean.fastq {output.barcodes} {output.ligation}) &> {log}
        # (pigz -f fastq/cleaned/{wildcards.sample}_tss-seq-clean.fastq) &>> {log}

rule fastqc_cleaned:
    input:
        "fastq/cleaned/{sample}_tss-seq-clean.fastq.gz"
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/cleaned/{sample}_clean_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/cleaned/fastqc_cleaned-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/cleaned) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/cleaned {input}) &>> {log}
        """

rule fastqc_aligned:
    input:
        lambda wc: "alignment/" + wc.sample + "_tss-seq-noPCRduplicates.bam" if wc.fqtype=="aligned_noPCRdup" else "alignment/" + wc.sample + "/unmapped.bam",
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/{fqtype}/{sample}_{fqtype}_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/{fqtype}/fastqc_{fqtype}-{sample}.log"
    wildcard_constraints:
        fqtype="aligned_noPCRdup|unaligned"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.fqtype}) &> {log}
        (bedtools bamtofastq -fq qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}_{wildcards.fqtype}.fastq -i {input}) &>> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}_{wildcards.fqtype}.fastq) &>> {log}
        (rm qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}_{wildcards.fqtype}.fastq) &>> {log}
        """

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}/{fname}/fastqc_data.txt", zip, sample=SAMPLES, fname=[os.path.split(v["fastq"])[1].split(".fastq")[0] + "_fastqc" for k,v in SAMPLES.items()]),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}_clean_fastqc/fastqc_data.txt", sample=SAMPLES),
        aligned_noPCRdup = expand("qual_ctrl/fastqc/aligned_noPCRdup/{sample}_aligned_noPCRdup_fastqc/fastqc_data.txt", sample=SAMPLES),
        unaligned = expand("qual_ctrl/fastqc/unaligned/{sample}_unaligned_fastqc/fastqc_data.txt", sample=SAMPLES),
    output:
        'qual_ctrl/fastqc/tss-seq-per_base_quality.tsv',
        'qual_ctrl/fastqc/tss-seq-per_tile_quality.tsv',
        'qual_ctrl/fastqc/tss-seq-per_sequence_quality.tsv',
        'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.tsv',
        'qual_ctrl/fastqc/tss-seq-per_sequence_gc.tsv',
        'qual_ctrl/fastqc/tss-seq-per_base_n.tsv',
        'qual_ctrl/fastqc/tss-seq-sequence_length_distribution.tsv',
        'qual_ctrl/fastqc/tss-seq-sequence_duplication_levels.tsv',
        'qual_ctrl/fastqc/tss-seq-adapter_content.tsv',
        'qual_ctrl/fastqc/tss-seq-kmer_content.tsv'
    run:
        shell("rm -f {output}")
        #for each statistic
        for outpath, stat, header in zip(output, ["Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Total Deduplicated Percentage", "Adapter Content", "Kmer Content"], ["base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus", "tile\tbase\tmean\tsample\tstatus",
        "quality\tcount\tsample\tstatus", "base\tg\ta\tt\tc\tsample\tstatus", "gc_content\tcount\tsample\tstatus", "base\tn_count\tsample\tstatus", "length\tcount\tsample\tstatus", "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus", "position\tpct\tsample\tstatus",
        "sequence\tcount\tpval\tobs_over_exp_max\tmax_position\tsample\tstatus" ]):
            for input_type in ["raw", "cleaned", "aligned_noPCRdup", "unaligned"]:
                for sample_id, fqc in zip(SAMPLES.keys(), input[input_type]):
                    shell("""awk -v sample_id={sample_id} -v input_type={input_type} 'BEGIN{{FS=OFS="\t"}} /{stat}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{print $0, sample_id, input_type}}' {fqc} | tail -n +2 >> {outpath}""")
            shell("""sed -i "1i {header}" {outpath}""")

rule plot_fastqc_summary:
    input:
        seq_len_dist = 'qual_ctrl/fastqc/tss-seq-sequence_length_distribution.tsv',
        per_tile = 'qual_ctrl/fastqc/tss-seq-per_tile_quality.tsv',
        per_base_qual = 'qual_ctrl/fastqc/tss-seq-per_base_quality.tsv',
        per_base_seq = 'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.tsv',
        per_base_n = 'qual_ctrl/fastqc/tss-seq-per_base_n.tsv',
        per_seq_gc = 'qual_ctrl/fastqc/tss-seq-per_sequence_gc.tsv',
        per_seq_qual = 'qual_ctrl/fastqc/tss-seq-per_sequence_quality.tsv',
        adapter_content = 'qual_ctrl/fastqc/tss-seq-adapter_content.tsv',
        seq_dup = 'qual_ctrl/fastqc/tss-seq-sequence_duplication_levels.tsv',
        kmer = 'qual_ctrl/fastqc/tss-seq-kmer_content.tsv'
    output:
        seq_len_dist = 'qual_ctrl/fastqc/tss-seq-sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/tss-seq-per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/tss-seq-per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/tss-seq-per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/tss-seq-per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/tss-seq-adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/tss-seq-sequence_duplication_levels.svg',
    script: "scripts/fastqc_summary.R"

#align to combined genome with Tophat2, WITHOUT reference transcriptome (i.e., the -G gff)
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
        fastq = "fastq/cleaned/{sample}_tss-seq-clean.fastq.gz"
    output:
        aligned = "alignment/{sample}/accepted_hits.bam",
        unaligned = "alignment/{sample}/unmapped.bam",
        summary = "alignment/{sample}/align_summary.txt",
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
    log: "logs/align/align_{sample}.log"
    shell:
        """
        (tophat2 --read-mismatches {params.read_mismatches} --read-gap-length {params.read_gap_length} --read-edit-dist {params.read_edit_dist} -o alignment/{wildcards.sample} --min-anchor-length {params.min_anchor_length} --splice-mismatches {params.splice_mismatches} --min-intron-length {params.min_intron_length} --max-intron-length {params.max_intron_length} --max-insertion-length {params.max_insertion_length} --max-deletion-length {params.max_deletion_length} --num-threads {threads} --max-multihits {params.max_multihits} --library-type fr-firststrand --segment-mismatches {params.segment_mismatches} --no-coverage-search --segment-length {params.segment_length} --min-coverage-intron {params.min_coverage_intron} --max-coverage-intron {params.max_coverage_intron} --min-segment-intron {params.min_segment_intron} --max-segment-intron {params.max_segment_intron} --b2-sensitive {params.idx_path}/{params.basename} {input.fastq}) &> {log}
        """

rule select_unique_mappers:
    input:
        "alignment/{sample}/accepted_hits.bam"
    output:
        temp("alignment/{sample}_tss-seq-uniquemappers.bam")
    threads: config["threads"]
    log: "logs/select_unique_mappers/select_unique_mappers_{sample}.log"
    shell: """
        (samtools view -buh -q 50 -@ {threads} {input} | samtools sort -T .{wildcards.sample} -@ {threads} - > {output}) &> {log}
        """

rule remove_PCR_duplicates:
    input:
        "alignment/{sample}_tss-seq-uniquemappers.bam"
    output:
        bam = "alignment/{sample}_tss-seq-noPCRduplicates.bam",
        bai = "alignment/{sample}_tss-seq-noPCRduplicates.bam.bai",
    log: "logs/remove_PCR_duplicates/remove_PCR_duplicates_{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output.bam}) &> {log}
        (samtools index {output.bam} > {output.bai}) &>> {log}
        """

rule read_processing_numbers:
    input:
        adapter = expand("logs/clean_reads/remove_adapter-{sample}.log", sample=SAMPLES),
        qual_trim = expand("logs/clean_reads/remove_3p_bc_and_trim-{sample}.log", sample=SAMPLES),
        align = expand("alignment/{sample}/align_summary.txt", sample=SAMPLES),
        nodups = expand("alignment/{sample}_noPCRdup.bam", sample=SAMPLES)
    output:
        "qual_ctrl/read_processing_summary.tsv"
    log: "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tadapter_removed\tquality_trimmed\tmapped\tunique_map\tnoPCRdup" > {output}) &> {log}""")
        for sample, adapter, qual_trim, align, nodups in zip(SAMPLES.keys(), input.adapter, input.qual_trim, input.align, input.nodups):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(grep -e "Reads written" {qual_trim} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(awk 'BEGIN{{ORS="\t"}} NR==3 || NR==4{{print $3}}' {align} >> {output}) &> {log}""")
            shell("""(samtools view -c {nodups} | awk '{{print $1}}' >> {output}) &> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$6=$5-$6; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing_summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing-loss.svg",
    script: "scripts/processing_summary.R"

rule genome_coverage:
    input:
        "alignment/{sample}-tss-seq-noPCRduplicates.bam"
    params:
        prefix = lambda wc: config["combinedgenome"]["experimental_prefix"] if wc.counttype=="counts" else config["combinedgenome"]["spikein_prefix"]
    output:
        # plmin = "coverage/{counttype}/{sample}_tss-seq-{counttype}-plmin.bedgraph",
        plus = "coverage/{counttype}/{sample}_tss-seq-{counttype}-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}_tss-seq-{counttype}-minus.bedgraph"
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/genome_coverage/genome_coverage-{sample}.log"
    shell: """
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | sed -n 's/{params.prefix}//p' | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | sed -n 's/{params.prefix}//p' | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        """
        # (genomeCoverageBed -bga -5 -ibam {input} | sed -n 's/{params.prefix}//p' | sort -k1,1 -k2,2n > {output.plmin}) &>> {log}

rule normalize:
    input:
        counts = "coverage/counts/{sample}-tss-counts-{strand}.bedgraph",
        plmin = lambda wc: "coverage/counts/" + wc.sample + "-tss-counts-plmin.bedgraph" if wc.norm=="libsizenorm" else "coverage/sicounts/" + wc.sample + "-tss-sicounts-plmin.bedgraph"
    params:
        scalefactor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
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
        plmin = expand("coverage/counts/{sample}-tss-counts-plmin.bedgraph", sample=sisamples),
        SIplmin = expand("coverage/sicounts/{sample}-tss-sicounts-plmin.bedgraph", sample=sisamples)
    params:
        group = [v["group"] for k,v in sisamples.items()]
    output:
        "qual_ctrl/all/spikein-counts.tsv"
    log: "logs/get_si_pct.log"
    run:
        shell("rm -f {output}")
        for name, exp, si, g in zip(SAMPLES.keys(), input.plmin, input.SIplmin, params.group):
            shell("""(echo -e "{name}\t{g}\t" $(awk 'BEGIN{{FS=OFS="\t"; ex=0; si=0}}{{if(NR==FNR){{si+=$4}} else{{ex+=$4}}}} END{{print ex+si, ex, si}}' {si} {exp}) >> {output}) &> {log}""")

#NOTE: since the conditiongroups_si and controlgroups_si lists are subsets of validgroups_si,
#in the case where there is only one passing sample, the group won't be included as is right now
rule plot_si_pct:
    input:
        "qual_ctrl/all/spikein-counts.tsv"
    output:
        plot = "qual_ctrl/{status}/{status}-spikein-plots.svg",
        stats = "qual_ctrl/{status}/{status}-spikein-stats.tsv"
    params:
        samplelist = lambda wc : list(sisamples.keys()) if wc.status=="all" else list(sipassing.keys()),
        conditions = conditiongroups_si,
        controls = controlgroups_si,
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
        lambda wc : FIGURES[wc.figure]["annotations"][wc.annotation]["path"]
    output:
        "{annopath}/stranded/{figure}_{annotation}-STRANDED.{ext}"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

def selectchrom(wc):
    if wc.strand in ["plus", "minus"]:
        if wc.norm=="sicounts":
            return config["genome"]["sichrsizes"]
        return config["genome"]["chrsizes"]
    if wc.norm=="sicounts":
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

rule compute_matrix:
    input:
        annotation = lambda wc: os.path.dirname(FIGURES[wc.figure]["annotations"][wc.annotation]["path"]) + "/stranded/" + wc.figure + "_" + wc.annotation + "-STRANDED" + os.path.splitext(FIGURES[wc.figure]["annotations"][wc.annotation]["path"])[1],
        bw = "coverage/{norm}/{sample}-tss-{norm}-{strand}.bw"
    output:
        dtfile = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}.tsv"),
        melted = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}-melted.tsv.gz")
    params:
        group = lambda wc : SAMPLES[wc.sample]["group"],
        refpoint = lambda wc: "TSS" if FIGURES[wc.figure]["parameters"]["type"]=="scaled" else FIGURES[wc.figure]["parameters"]["refpoint"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        binsize = lambda wc: FIGURES[wc.figure]["parameters"]["binsize"],
        binstat = lambda wc: FIGURES[wc.figure]["parameters"]["binstat"],
        nan_afterend = lambda wc: [] if FIGURES[wc.figure]["parameters"]["type"]=="scaled" or not FIGURES[wc.figure]["parameters"]["nan_afterend"] else "--nanAfterEnd",
        anno_label = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["label"]
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{annotation}_{sample}_{norm}-{strand}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} -g {params.group} -s {wildcards.sample} -a {params.anno_label} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}-melted.tsv.gz", annotation=[k for k,v in FIGURES[wc.figure]["annotations"].items()], sample=SAMPLES, figure=wc.figure, norm=wc.norm, strand=wc.strand)
    output:
        "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{norm}-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{figure}_{norm}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        matrices = expand("datavis/{{figure}}/{{norm}}/{{figure}}-allsamples-allannotations-{{norm}}-{strand}.tsv.gz", strand=["SENSE", "ANTISENSE"]),
        annotations = lambda wc: [v["path"] for k,v in FIGURES[wc.figure]["annotations"].items()]
    output:
        heatmap_sample_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample-bothstrands.svg",
        heatmap_sample_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample-sense.svg",
        heatmap_sample_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample-antisense.svg",
        heatmap_group_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup-bothstrands.svg",
        heatmap_group_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup-sense.svg",
        heatmap_group_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup-antisense.svg",
        metagene_sample_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-bothstrands.svg",
        metagene_sample_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-sense.svg",
        metagene_sample_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-antisense.svg",
        metagene_sample_overlay_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-bothstrands.svg",
        metagene_sample_overlay_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-sense.svg",
        metagene_sample_overlay_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-antisense.svg",
        metagene_group_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-bothstrands.svg",
        metagene_group_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-sense.svg",
        metagene_group_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-antisense.svg",
        metagene_sampleclust_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-bothstrands.svg",
        metagene_sampleclust_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-sense.svg",
        metagene_sampleclust_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-antisense.svg",
        metagene_groupclust_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-bothstrands.svg",
        metagene_groupclust_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-sense.svg",
        metagene_groupclust_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-antisense.svg",
    params:
        # abusing snakemake a bit here...using params as output paths in order to use lambda functions
        annotations_out = lambda wc: ["datavis/" + wc.figure + "/" + wc.norm + "/" + wc.condition + "-v-" + wc.control + "/" + wc.status + "/" + annotation + "_cluster-" + str(cluster) + ".bed" for annotation in FIGURES[wc.figure]["annotations"] for cluster in range(1, FIGURES[wc.figure]["annotations"][annotation]["n_clusters"]+1)],
        clusters_out = lambda wc: ["datavis/" + wc.figure + "/" + wc.norm + "/" + wc.condition + "-v-" + wc.control + "/" + wc.status + "/" + annotation + ".pdf" for annotation in FIGURES[wc.figure]["annotations"]],
        samplelist = plotcorrsamples,
        plottype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        pct_cutoff = lambda wc: FIGURES[wc.figure]["parameters"]["pct_cutoff"],
        log_transform = lambda wc: str(FIGURES[wc.figure]["parameters"]["log_transform"]).upper(),
        pcount = lambda wc: 0 if not FIGURES[wc.figure]["parameters"]["log_transform"] else FIGURES[wc.figure]["parameters"]["pseudocount"],
        trim_pct = lambda wc: FIGURES[wc.figure]["parameters"]["trim_pct"],
        refpointlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
        endlabel = lambda wc:  "HAIL SATAN" if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["three_prime_label"],
        cmap = lambda wc: FIGURES[wc.figure]["parameters"]["heatmap_colormap"],
        sortmethod = lambda wc: FIGURES[wc.figure]["parameters"]["arrange"],
        cluster_scale = lambda wc: "FALSE" if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else str(FIGURES[wc.figure]["parameters"]["cluster_scale"]).upper(),
        cluster_samples = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else cluster_samples(wc.status, wc.norm, FIGURES[wc.figure]["parameters"]["cluster_conditions"], FIGURES[wc.figure]["parameters"]["cluster_strands"]),
        cluster_five = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_five"],
        cluster_three = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_three"],
        k = lambda wc: [v["n_clusters"] for k,v in FIGURES[wc.figure]["annotations"].items()],
    script:
        "scripts/plot_tss_figures.R"

rule union_bedgraph:
    input:
        exp = expand("coverage/{{norm}}/{sample}-tss-{{norm}}-SENSE.bedgraph", sample=SAMPLES)
    output:
        exp = "coverage/{norm}/union-bedgraph-allsamples-{norm}.tsv.gz",
    params:
        names = list(SAMPLES.keys())
    log: "logs/union_bedgraph-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz > {output.exp}) &> {log}
        """

#TODO: plot correlations over multiple binsizes, as well as over annotations
rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-allsamples-{norm}.tsv.gz"
    output:
        "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-tss-{status}-{norm}-correlations.svg"
    params:
        pcount = 0.1,
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule call_tss_peaks:
    input:
        bw = lambda wc: "coverage/counts/" + wc.sample + "-tss-counts-SENSE.bw" if wc.type=="exp" else "coverage/sicounts/" + wc.sample + "-tss-sicounts-SENSE.bw"
    output:
        smoothed = expand("peakcalling/sample_peaks/{{sample}}-{{type}}-smoothed-bw{bandwidth}-{strand}.bw", strand=["plus","minus"], bandwidth = config["peakcalling"]["bandwidth"]),
        peaks = "peakcalling/sample_peaks/{sample}-{type}-allpeaks.narrowPeak"
    params:
        name = lambda wc: wc.sample + "-" + wc.type,
        bandwidth = config["peakcalling"]["bandwidth"],
        window = config["peakcalling"]["local-bg-window"]
    log: "logs/call_tss_peaks/call_tss_peaks-{sample}-{type}.log"
    shell: """
        (python scripts/tss-peakcalling.py -i {input.bw} -n {params.name} -w {params.window} -b {params.bandwidth} -o peakcalling/sample_peaks) &> {log}
        """

rule tss_peaks_idr:
    input:
        #NOTE: for now we take the first two samples since the IDR script only takes two
        #change this if we find a better way to aggregate results
        lambda wc: ["peakcalling/sample_peaks/" + x + "-" + wc.type + "-allpeaks.narrowPeak" for x in PASSING if PASSING[x]['group']==wc.group][0:2]
    output:
        allpeaks = "peakcalling/{group}/{group}-{type}-idrpeaks-all.tsv",
        filtered = "peakcalling/{group}/{group}-{type}-idrpeaks-filtered.tsv",
        narrowpeak = "peakcalling/{group}/{group}-{type}-idrpeaks.narrowPeak",
        summits = "peakcalling/{group}/{group}-{type}-idrpeaks-summits.bed",
    params:
        idr = int(-125*log2(config["peakcalling"]["idr"]))
    log: "logs/tss_peaks_idr/tss_peaks_idr-{group}-{type}.log"
    shell: """
        idr -s {input} --input-file-type narrowPeak --rank q.value -o {output.allpeaks} -l {log} --plot --peak-merge-method max
        (awk '$5>{params.idr} || $9=="inf"' {output.allpeaks} | tee {output.filtered} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $11, $12, $10}}' | sed "s/-minus//g;s/-plus//g" | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.summits}) &> {log}
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
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-genic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | cut --complement -f11-13,15-17 > {output}
        """

rule classify_peaks_intragenic:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-intragenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}} $6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule classify_peaks_antisense:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        transcripts = config["genome"]["transcripts"]
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-antisense.tsv"
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
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-convergent.tsv"
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
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-divergent.tsv"
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
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-intergenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.transcripts} {input.orfs} -wa -v | bedtools intersect -a stdin -b {input.genic_anno} -wa -v -s | bedtools intersect -a stdin -b {input.annotation} -wa > {output}
        """

rule peakstats:
    input:
        expand("peakcalling/{group}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
    output:
        table = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaknumbers.tsv",
        size = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaksizes-histogram.svg",
        violin_area = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaksizes-violin-equalarea.svg",
        violin_count = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaksizes-violin-countscaled.svg",
        dist = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peakdistances.svg"
    params:
        groups = lambda wc: [g for sublist in zip(controlgroups, conditiongroups) for g in sublist] if wc.condition=="all" else [wc.control, wc.condition]
    script:
        "scripts/peakstats.R"

rule combine_tss_peaks:
    input:
        cond = "peakcalling/{condition}/{condition}-{type}-idrpeaks-filtered.tsv",
        ctrl = "peakcalling/{control}/{control}-{type}-idrpeaks-filtered.tsv",
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-peaks.bed"
    shell: """
        sort -k1,1 -k2,2n {input.cond} | bedtools multiinter -i stdin <(sort -k1,1 -k2,2n {input.ctrl}) -cluster | cut -f1-3 | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule map_counts_to_peaks:
    input:
        bed = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-peaks.bed",
        bg = lambda wc: "coverage/counts/" + wc.sample + "-tss-counts-SENSE.bedgraph" if wc.type=="exp" else "coverage/sicounts/" + wc.sample + "-tss-sicounts-SENSE.bedgraph"
    output:
        temp("diff_exp/{condition}-v-{control}/{sample}_tss-seq-{type}-peakcounts.tsv")
    log: "logs/map_counts_to_peaks/map_counts_to_peaks-{condition}-v-{control}_{sample}-{type}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-"$2"-"$3, $4}}' &> {output}) &> {log}
        """

def getsamples(ctrl, cond):
    return [k for k,v in PASSING.items() if v["group"] in [ctrl, cond]]

rule combine_peak_counts:
    input:
        lambda wc : ["diff_exp/" + wc.condition + "-v-" + wc.control + "/" + x + "_tss-seq-" + wc.type + "-peakcounts.tsv" for x in getsamples(wc.control, wc.condition)]
    output:
        "diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-{type}-peak-counts.tsv"
    params:
        n = lambda wc: 2*len(getsamples(wc.control, wc.condition)),
        names = lambda wc: "\t".join(getsamples(wc.control, wc.condition))
    log: "logs/get_peak_counts/get_peak_counts-{condition}-v-{control}_{type}.log"
    shell: """
        (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
        """

rule differential_expression:
    input:
        expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-experimental-peak-counts.tsv",
        sicounts = lambda wc: "diff_exp/" + wc.condition + "-v-" + wc.control + "/" + wc.condition + "-v-" + wc.control + "_tss-seq-spikein-peak-counts.tsv" if wc.norm=="spikenorm" else []
    params:
        samples = lambda wc : getsamples(wc.control, wc.condition),
        groups = lambda wc : [PASSING[x]["group"] for x in getsamples(wc.control, wc.condition)],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"])
    output:
        results_all = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-all.tsv",
        results_up = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-up.tsv",
        results_down = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-down.tsv",
        results_unch = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-unchanged.tsv",
        normcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-counts-sizefactornorm.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-counts-rlogtransformed.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-qc-plots.svg"
    script:
        "scripts/call_de_peaks.R"

rule diffexp_results_to_narrowpeak:
    input:
        condition_coverage = lambda wc: expand("coverage/{norm}/{sample}-tss-{norm}-SENSE.bw", sample=getsamples(wc.condition, wc.condition), norm=wc.norm),
        control_coverage = lambda wc: expand("coverage/{norm}/{sample}-tss-{norm}-SENSE.bw", sample=getsamples(wc.control, wc.control), norm=wc.norm),
        diffexp_results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        summit_bed = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}-summits.bed",
    log: "logs/diffexp_results_to_narrowpeak/diffexp_results_to_narrowpeak-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/diffexp_results_to_narrowpeak.py -i {input.condition_coverage} -j {input.control_coverage} -d {input.diffexp_results} -n {output.narrowpeak} -b {output.summit_bed}) &> {log}
        """

rule classify_genic_diffexp_peaks:
    input:
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}-summits.bed",
    log : "logs/classify_genic_diffexp_peaks/classify_genic_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b <(cut -f1-6 {input.annotation}) -wo -s | cut --complement -f21 | cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_score\ttranscript_name\ttranscript_score\ttranscript_strand")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.annotation} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """
    # shell: """
    #     (bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10}} $6=="-"{{print $4, $8, $9, $10}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name") - | tee {output.table} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $2, $4, $5, $1, $7":"$8, $3}}' > {output.bed}) &> {log}
    #     """

rule get_de_intragenic:
    input:
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        orf_anno = config["genome"]["orf-annotation"],
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.tsv"
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic-summits.bed"
    log: "logs/get_de_intragenic/get_de_intragenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.orfs}) -wo -s | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=summit-$12}} $6=="-"{{$17=$13-(summit+1)}} {{print $0}}'  | cat <(paste <(head -n 1 {input.results}) <(echo -e "orf_chrom\torf_start\torf_end\torf_score\torf_name\torf_score\torf_strand\tatg_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """
    # shell: """
    #     (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}} $6=="-"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\tORF_start\tORF_end\tORF_name\tdist_ATG_to_peak") - | tee {output.table} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $2, $4, $5, $1, $7":"$8, $3}}' > {output.bed}) &> {log}
    #     """

rule get_de_intragenic_frequency:
    input:
        orfs = config["genome"]["orf-annotation"],
        intrabed = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic.bed"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
    shell: """
        bedtools intersect -a {input.orfs} -b {input.intrabed} -c -s > {output}
        """

rule plot_de_intragenic_frequency:
    input:
        "diff_exp/{condition}-v-{control}/{norm}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
    output:
        "diff_exp/{condition}-v-{control}/{norm}/intragenic/intrafreq/{condition}-v-{control}-intragenic-{norm}-{direction}-freqperORF.svg"
    log: "logs/plot_intragenic_frequency/plot_intragenic_frequency-{condition}-v-{control}-{norm}-{direction}.log"
    script: "scripts/intrafreq.R"

rule get_de_antisense:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        transcripts = config["genome"]["transcripts"],
        totalresults = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.tsv"
    output:
        table = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}-results-{norm}-{direction}-antisense.tsv",
        bed = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}-results-{norm}-{direction}-antisense.bed"
    log : "logs/get_de_antisense/get_de_antisense-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcripts} -wo -S | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}} $6=="-"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name\tdist_peak_to_senseTSS") - | tee {output.table} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $2, $4, $5, $1, $7":"$8, $3}}' > {output.bed}) &> {log}
        """

rule get_de_convergent:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.tsv"
    output:
        table = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}-results-{norm}-{direction}-convergent.tsv",
        bed = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}-results-{norm}-{direction}-convergent.bed"
    log : "logs/get_de_convergent/get_de_convergent-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}} $6=="-"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name\tdist_peak_to_senseTSS") - | tee {output.table} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $2, $4, $5, $1, $7":"$8, $3}}' > {output.bed}) &> {log}
        """

rule get_de_divergent:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.tsv"
    output:
        table = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}-results-{norm}-{direction}-divergent.tsv",
        bed = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}-results-{norm}-{direction}-divergent.bed"
    log : "logs/get_de_divergent/get_de_divergent-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}} $6=="-"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name\tdist_peak_to_senseTSS") - | tee {output.table} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $2, $4, $5, $1, $7":"$8, $3}}' > {output.bed}) &> {log}
        """


rule get_de_intergenic:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-{direction}.tsv"
    output:
        table = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}-results-{norm}-{direction}-intergenic.tsv",
        bed = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}-results-{norm}-{direction}-intergenic.bed",
    log : "logs/get_de_intergenic/get_de_intergenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcripts} {input.orfs} -wa -v | bedtools intersect -a stdin -b {input.genic_anno} -wa -v -s | bedtools intersect -a stdin -b {input.annotation} -wo | awk 'BEGIN{{FS=OFS="\t"}}$6=="+"{{print $4, $8, $9}}$6=="-"{{print $4, $8, $9}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\tregion_start\tregion_end") - | tee {output.table} | awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print $2, $4, $5, $1, $7":"$8, $3}}' > {output.bed}) &> {log}
        """

rule get_intra_orfs:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-de-clusters-{norm}-{direction}-intragenic.tsv",
        fasta = config["genome"]["fasta"]
    output:
        "diff_exp/{condition}-v-{control}/{norm}/intragenic/intragenic-orfs/{condition}-v-{control}-{norm}-{direction}-intragenic-orfs.tsv"
    params:
        max_upstr_atgs = config["max-upstr-atgs"],
        max_search_dist = 2000
    log: "logs/get_intra_orfs/get_intra_orfs-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (python scripts/find_intra_orfs.py -p {input.peaks} -f {input.fasta} -m {params.max_search_dist} -a {params.max_upstr_atgs} -o {output}) &> {log}
        """

#TODO: account for double-counted peaks when a peak overlaps more than one annotation (more than one genic region, for example)
rule summarise_de_results:
    input:
        total = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-results-{norm}-all.tsv",
        genic = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}-results-{norm}-all-genic.tsv",
        intragenic = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-all-intragenic.tsv",
        antisense = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}-results-{norm}-all-antisense.tsv",
        convergent = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}-results-{norm}-all-convergent.tsv",
        divergent = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}-results-{norm}-all-divergent.tsv",
        intergenic = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}-results-{norm}-all-intergenic.tsv",
    params:
        lfc = config["deseq"]["fold-change-threshold"],
        alpha = config["deseq"]["fdr"]
    output:
        summary_table = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-summary.tsv",
        summary = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-summary.svg",
        maplot = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-maplot.svg",
        volcano = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-volcano.svg",
        volcano_free = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}-{norm}-diffexp-volcano-freescale.svg",
    script: "scripts/de_summary.R"

rule genic_v_class:
    input:
        genic = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}-results-{norm}-all-genic.tsv",
        intragenic = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-all-intragenic.tsv",
        antisense = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}-results-{norm}-all-antisense.tsv",
        convergent = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}-results-{norm}-all-convergent.tsv",
        divergent = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}-results-{norm}-all-divergent.tsv",
    params:
        path = "diff_exp/{condition}-v-{control}/{norm}/genic_v_class/"
    output:
        figure = "diff_exp/{condition}-v-{control}/{norm}/genic_v_class/{condition}-v-{control}-{norm}-genic-v-class.svg",
        tables = expand("diff_exp/{{condition}}-v-{{control}}/{{norm}}/genic_v_class/{{condition}}-v-{{control}}-{{norm}}-genic-v-{ttype}.tsv", ttype=["intragenic", "antisense", "convergent", "divergent"])
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
#0. filter out double counted peaks (sometimes a peak can be 'genic' for two genes, causing it to be listed twice)
#1. with the START of the peak as reference, extend annotation to upstream and 'downstream' distances
#2. if multiple annotations overlap on same strand, keep the one that is the most significant (avoid multiple-counting, so that poorly called peaks which are erroneously split into multiple peaks do not bias the motif enrichment analyses.)
#3. intersect with motif file
rule get_upstream_motifs:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed",
        chrsizes = config["genome"]["chrsizes"],
        motifs = "motifs/allmotifs.bed"
    output:
        "motifs/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-{direction}-{category}-motifs.tsv.gz"
    params:
        upstr = config["motifs"]["enrichment-upstream"],
        dnstr = config["motifs"]["enrichment-downstream"]
    log: "logs/get_upstream_motifs/get_upstream_motifs-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (uniq {input.peaks} | bedtools flank -l {params.upstr} -r 0 -s -i stdin -g {input.chrsizes} | bedtools slop -l 0 -r {params.dnstr} -s -i stdin -g {input.chrsizes} | bedtools cluster -s -d 0 -i stdin | sed 's/:/\t/g' | bedtools groupby -g 8 -c 6 -o max -full -i stdin | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5":"$6, $7}}' | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b {input.motifs} -sorted -wao | awk 'BEGIN{{FS="\t|:"; OFS="\t"}}{{print $1, $4, $5, $6, $11, $14, $9, $10, $12}}' | cat <(echo -e "chrom\ttss_peak_id\tpeak_lfc\tpeak_logpadj\tmotif_id\tmotif_alt_id\tmotif_start\tmotif_end\tmotif_logpadj") - | pigz -f > {output}) &> {log}
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
    log: "logs/get_random_motifs/get_random_motifs.log"
    shell: """
        (bedtools random -l {params.size} -n {params.n} -g {input.chrsizes} | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b {input.motifs} -sorted -wao | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $4, 0, 0, $10, $13, $8, $9, $11}}' | cat <(echo -e "chrom\ttss_peak_id\tpeak_lfc\tpeak_logpadj\tmotif_id\tmotif_alt_id\tmotif_start\tmotif_end\tmotif_logpadj") - | pigz -f > {output}) &> {log}
        """

rule test_motif_enrichment:
    input:
        fimo_pos = "motifs/{condition}-v-{control}/{norm}/{condition}-v-{control}_{norm}-{direction}-{category}-motifs.tsv.gz",
        fimo_neg = lambda wc: "motifs/" + wc.condition + "-v-" + wc.control + "/" + wc.norm + "/" + wc.condition + "-v-" + wc.control + "_" + wc.norm + "-unchanged-" + wc.category + "-motifs.tsv.gz" if wc.negative=="unchanged" else "motifs/random-motifs.tsv.gz"
    params:
        pval_cutoff = config["motifs"]["fimo-pval"],
        alpha= config["motifs"]["enrichment-fdr"],
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated"
    output:
        tsv = "motifs/{condition}-v-{control}/{norm}/{negative}/{condition}-v-{control}_{norm}-{direction}-v-{negative}-{category}-motif_enrichment.tsv",
        plot = "motifs/{condition}-v-{control}/{norm}/{negative}/{condition}-v-{control}_{norm}-{direction}-v-{negative}-{category}-motif_enrichment.svg",
    script: "scripts/motif_enrichment.R"

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

#TODO: fix the statistical test
rule peak_positioning:
    input:
        peaks = "diff_exp/{condition}-v-{control}/{norm}/{ttype}/{condition}-v-{control}-results-{norm}-{direction}-{ttype}.tsv",
    params:
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated",
        n_bins = 100
    output:
        relative = "diff_exp/{condition}-v-{control}/{norm}/{ttype}/{condition}-v-{control}-relative-distances-{norm}-{direction}-{ttype}.svg",
        absolute = "diff_exp/{condition}-v-{control}/{norm}/{ttype}/{condition}-v-{control}-absolute-distances-{norm}-{direction}-{ttype}.svg",
        fc_signif = "diff_exp/{condition}-v-{control}/{norm}/{ttype}/{condition}-v-{control}-foldchange-significance-{norm}-{direction}-{ttype}.svg"
    script:
        "scripts/length_bias.R"

#TODO: move sequence analysis to separate pipeline
rule get_gc_percentage:
    input:
        fasta = config["genome"]["fasta"],
    params:
        binsize = 11 #must be odd integer
    output:
        os.path.splitext(config["genome"]["fasta"])[0] + "-GC_pct.bw"
    shell: """
        python scripts/gc_content.py -f {input.fasta} -w {params.binsize} -o {output}
        """

rule gene_ontology:
    input:
        universe = lambda wc: config["genome"]["orf-annotation"] if wc.category=="intragenic" else config["genome"]["transcripts"],
        diffexp_path = "diff_exp/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.tsv",
        go_anno_path = config["gene_ontology_mapping_file"]
    params:
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated" if wc.direction=="down" else wc.direction
    output:
        results = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-GO-{norm}-{direction}-{category}-results.tsv",
        enriched_combined = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-GO-{norm}-{direction}-{category}-enriched-all.svg",
        depleted_combined = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-GO-{norm}-{direction}-{category}-depleted-all.svg",
        enriched_facet = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-GO-{norm}-{direction}-{category}-enriched-facetted.svg",
        depleted_facet = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}-GO-{norm}-{direction}-{category}-depleted-facetted.svg",
    script:
        "scripts/gene_ontology.R"

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

rule get_seqlogo_data:
    input:
        bam = lambda wc: expand("alignment/{sample}-noPCRdup.bam", sample=getsamples(wc.group, wc.group)),
        bed = lambda wc: "peakcalling/" + wc.category + "/" + wc.group + "-exp-idrpeaks-" + wc.category + ".tsv" if wc.category != "all" else [],
        chrsizes = config["genome"]["chrsizes"],
        fasta = config["genome"]["fasta"]
    params:
        prefix = config["combinedgenome"]["experimental_prefix"],
        slop = int(config["consensus"]["window"]),
        windowsize = int(2*config["consensus"]["window"]+1),
        intersect_cmd = lambda wc: "bedtools intersect -wa -s -a stdin -b peakcalling/" + wc.category + "/" + wc.group + "-exp-idrpeaks-" + wc.category + ".tsv | " if wc.category != "all" else ""
    output:
        seqlogo_data = "seq_logos/{group}/{group}-{category}-seqlogo.tsv",
    log: "logs/get_seqlogo_data/get_seqlogo_data-{group}-{category}.log"
    run:
        fasta_stats = subprocess.run(args="fasta-get-markov " + input.fasta + " | tail -n +4", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()
        composition = "{" + ",".join(["'" + fasta_stats[x] + "':" + str(1e2*float(fasta_stats[x+1])) for x in range(0,8,2)]) + "}"
        shell("""(samtools merge - {input.bam} | bedtools bamtobed -i stdin | grep -e "{params.prefix}" | sed 's/{params.prefix}//g' | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6}}' | {params.intersect_cmd} bedtools slop -b {params.slop} -s -i stdin -g {input.chrsizes} | awk '$3-$2=={params.windowsize}' | bedtools getfasta -s -fi {input.fasta} -bed stdin | weblogo --format logodata --sequence-type dna --composition "{composition}" > {output.seqlogo_data}) &> {log}""")

rule plot_seqlogos:
    input:
        seqlogo_data = expand("seq_logos/{{group}}/{{group}}-{category}-seqlogo.tsv", category = ["all"] + CATEGORIES),
        fasta = config["genome"]["fasta"]
    params:
        tss_classes = ["all"] + CATEGORIES,
        slop = int(config["consensus"]["window"]),
    output:
        "seq_logos/{group}/{group}-seqlogos.svg",
    log: "logs/plot_seqlogos/plot_seqlogos-{group}.log"
    run:
        gc_pct = float(subprocess.run(args="fasta-get-markov " + input.fasta + " | tail -n +4", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()[3])*2
        shell(""" (Rscript scripts/plot_seqlogo.R --input {input.seqlogo_data} --tss_classes {params.tss_classes} -b {params.slop} --gc_pct {gc_pct} -l {wildcards.group} -o {output}) &> {log}""")

