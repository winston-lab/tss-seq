#!/usr/bin/env python
import os
import re
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
    build_library_size_table, plot_si_pct,
    make_stranded_genome, make_stranded_annotations,
    build_genic_annotation, build_convergent_annotation, build_divergent_annotation, build_intergenic_annotation,
    classify_peaks_genic, classify_peaks_intragenic, classify_peaks_antisense,
    classify_peaks_convergent, classify_peaks_divergent, classify_peaks_intergenic,
    combine_tss_peaks, map_counts_to_peaks, combine_peak_counts,
    classify_genic_diffexp_peaks, classify_intragenic_diffexp_peaks, classify_antisense_diffexp_peaks,
    classify_convergent_diffexp_peaks, classify_divergent_diffexp_peaks, classify_intergenic_diffexp_peaks,
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
        'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=SAMPLES),
        ##quality controls
        #"qual_ctrl/read_processing-loss.svg",
        #expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all", "passing"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-tss-{{status}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status = ["all", "passing"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-tss-{{status}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status = ["all", "passing"]),
        #coverage
        expand("coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bw", norm=["spikenorm","libsizenorm", "counts", "sicounts"], sample=SAMPLES, strand=["SENSE","ANTISENSE","plus","minus"]),
        #peakcalling on all samples
        expand("peakcalling/sample_peaks/{sample}_experimental-allpeaks.narrowPeak", sample=SAMPLES),
        expand("peakcalling/sample_peaks/{sample}_spikein-allpeaks.narrowPeak", sample=sisamples),
        ##IDR for all groups which have at least two passing samples
        expand("peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak", group=validgroups),
        expand("peakcalling/{group}/{group}_spikein-idrpeaks.narrowPeak", group=validgroups_si),
        ##classify peaks into categories
        #expand("peakcalling/{group}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
        #expand("peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peakdistances.svg", zip, condition=conditiongroups + ["all"], control=controlgroups + ["all"]),
        #combine called peaks for conditions vs control
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_experimental-peaks.bed", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_spikein-peaks.bed", zip, condition=conditiongroups_si, control=controlgroups_si),
        #differential expression of peaks
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-experimental-peak-counts.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-spikein-peak-counts.tsv", zip, condition=conditiongroups_si, control=controlgroups_si),
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-results-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups),direction=["all", "up", "unchanged", "down"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{condition}-v-{control}_tss-seq-spikenorm-diffexp-results-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si),direction=["all", "up", "unchanged", "down"]),
        #categorize differentially expressed peaks
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-results-{{category}}-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups), direction = ["all","up","unchanged","down"], category=CATEGORIES),
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

rule build_library_size_table:
    input:
        bams = expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=sisamples),
        bais = expand("alignment/{sample}_tss-seq-noPCRduplicates.bam.bai", sample=sisamples),
        chrsizes = config["combinedgenome"]["chrsizes"]
    output:
        "qual_ctrl/all/tss-seq-library-sizes.tsv"
    params:
        groups = [v["group"] for k,v in sisamples.items()],
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
        spikein_prefix = config["combinedgenome"]["spikein_prefix"],
    log: "logs/build_library_size_table.log"
    run:
        shell("""(echo -e "sample\tgroup\ttotal_counts\texperimental_counts\tspikein_counts" > {output}) &> {log} """)
        for sample, group, bam in zip(sisamples.keys(), params.groups, input.bams):
            shell("""(paste <(echo -e "{sample}\t{group}") <(samtools view -c {bam}) <(grep -oh "\w*{params.exp_prefix}\w*" {input.chrsizes} | xargs samtools view -c {bam}) <(samtools view -c {bam}) <(grep -oh "\w*{params.spikein_prefix}\w*" {input.chrsizes} | xargs samtools view -c {bam}) >> {output}) &>> {log}""")

#NOTE: since the conditiongroups_si and controlgroups_si lists are subsets of validgroups_si,
#in the case where there is only one passing sample, the group won't be included as is right now
#TODO: fix condition-v-control comparisons
rule plot_si_pct:
    input:
        "qual_ctrl/all/tss-seq-library-sizes.tsv"
    output:
        plot = "qual_ctrl/{status}/tss-seq-spikein-plots-{status}.svg",
        stats = "qual_ctrl/{status}/tss-seq-spikein-stats-{status}.tsv"
    params:
        samplelist = lambda wc : list(sisamples.keys()) if wc.status=="all" else list(sipassing.keys()),
        conditions = conditiongroups_si,
        controls = controlgroups_si,
    script: "scripts/plotsipct.R"

#make 'stranded' genome for datavis purposes
rule make_stranded_genome:
    input:
        exp = config["genome"]["chrsizes"],
        spikein = config["genome"]["sichrsizes"]
    output:
        exp = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
        spikein = os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv",
    log: "logs/make_stranded_genome.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.exp} > {output.exp}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.spikein} > {output.spikein}) &>> {log}
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

def getsamples(ctrl, cond):
    return [k for k,v in PASSING.items() if v["group"] in [ctrl, cond]]

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

include: "rules/tss-seq_clean_reads.smk"
include: "rules/tss-seq_alignment.smk"
include: "rules/tss-seq_fastqc.smk"
include: "rules/tss-seq_genome_coverage.smk"
include: "rules/tss-seq_datavis.smk"
include: "rules/tss-seq_peakcalling.smk"
include: "rules/tss-seq_differential_expression.smk"
include: "rules/tss-seq_build_annotations.smk"
include: "rules/tss-seq_classify_peaks.smk"
include: "rules/tss-seq_gene_ontology.smk"
include: "rules/tss-seq_motifs.smk"
include: "rules/tss-seq_sequence_logos.smk"

