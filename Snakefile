#!/usr/bin/env python

import os
import re
import subprocess
import itertools
from math import log2, log10

configfile: "config.yaml"

subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}

#groups which have at least two passing samples, so that they are valid for peakcalling and diff exp
validgroups = set([z for z in [PASSING[x]['group'] for x in PASSING] if [PASSING[x]['group'] for x in PASSING].count(z)>=2])
validgroups_si = set([z for z in [PASSING[x]['group'] for x in PASSING if PASSING[x]['spikein']=="y"] if [PASSING[x]['group'] for x in PASSING].count(z)>=2])
controlgroups = [v for k,v in config["comparisons"]["libsizenorm"].items() if v in validgroups]
conditiongroups = [k for k,v in config["comparisons"]["libsizenorm"].items() if k in validgroups]
controlgroups_si = [v for k,v in config["comparisons"]["spikenorm"].items() if v in validgroups_si]
conditiongroups_si = [k for k,v in config["comparisons"]["spikenorm"].items() if k in validgroups_si]

CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

FIGURES = config["figures"]

#get all motif names from motif databases, cleaning nasty characters in some motif names
MOTIFS = set(subprocess.run(args="meme2meme " + " ".join(config["motifs"]["databases"]) + " | grep -e '^MOTIF' | cut -d ' ' -f2 | sed 's/\//_/g; s/&/_/g; s/{/[/g; s/}/]/g' ", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split())

status_norm_sample_dict = {
    "all":
        {   "libsizenorm" : SAMPLES,
            "spikenorm" : SISAMPLES
        },
    "passing":
        {   "libsizenorm" : PASSING,
            "spikenorm" : SIPASSING
        }
    }

def get_samples(status, norm, groups):
    if "all" in groups:
        return(list(status_norm_sample_dict[status][norm].keys()))
    else:
        return([k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in groups])

def cluster_samples(status, norm, cluster_groups, cluster_strands):
    ll = []
    for group, strand in zip(cluster_groups, cluster_strands):
        sublist = [k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in cluster_groups]
        if strand in ["sense", "both"]:
            ll.append([f"{sample}-sense" for sample in sublist])
        if strand in ["antisense", "both"]:
            ll.append([f"{sample}-antisense" for sample in sublist])
    return(list(itertools.chain(*ll)))

include: "rules/tss-seq_clean_reads.smk"
include: "rules/tss-seq_alignment.smk"
include: "rules/tss-seq_fastqc.smk"
include: "rules/tss-seq_library-processing-summary.smk"
include: "rules/tss-seq_sample_similarity.smk"
include: "rules/tss-seq_genome_coverage.smk"
include: "rules/tss-seq_datavis.smk"
include: "rules/tss-seq_peakcalling.smk"
include: "rules/tss-seq_differential_expression.smk"
include: "rules/tss-seq_classify_peaks.smk"
include: "rules/tss-seq_gene_ontology.smk"
include: "rules/tss-seq_motifs.smk"
include: "rules/tss-seq_sequence_logos.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules:
    all,
    make_stranded_genome,
    genic_v_class,

rule all:
    input:
        #FastQC
        'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bw", norm=["spikenorm","libsizenorm", "counts", "sicounts"], sample=SAMPLES, strand=["SENSE","ANTISENSE","plus","minus"]),
        #quality controls
        "qual_ctrl/read_processing/tss-seq_read-processing-loss.svg",
        expand("qual_ctrl/spikein/tss-seq_spikein-plots-{status}.svg", status=["all", "passing"]),
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_tss-seq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all", "passing"], windowsize = config["scatterplot_binsizes"]),
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_tss-seq-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all", "passing"], windowsize = config["scatterplot_binsizes"]),
        #peakcalling on all samples
        expand("peakcalling/sample_peaks/{sample}_experimental-allpeaks.narrowPeak", sample=SAMPLES),
        expand("peakcalling/sample_peaks/{sample}_spikein-allpeaks.narrowPeak", sample=SISAMPLES),
        #IDR for all groups which have at least two passing samples
        expand("peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak", group=validgroups),
        expand("peakcalling/{group}/{group}_spikein-idrpeaks.narrowPeak", group=validgroups_si),
        #classify peaks into categories
        expand("peakcalling/{group}/{category}/{group}-experimental-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
        #some peak statistics (size, distance from sense TSS, ATG)
        expand("peakcalling/{group}/{group}-experimental-peak-stats.tsv", group=validgroups),
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
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_tss-seq-spikenorm-diffexp-results-{{category}}-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["all","up","unchanged","down"], category=CATEGORIES),
        #differential expression summary
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups, control=controlgroups), plot = ["mosaic", "maplot", "volcano"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{condition}-v-{control}_tss-seq-spikenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), plot = ["mosaic", "maplot", "volcano"]),
        #random distributions for differentially expressed peaks
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/intragenic/position_bias/{condition}-v-{control}_tss-seq-libsizenorm-intragenic-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down"], reference=["ATG", "RELATIVE", "STOP"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/intragenic/position_bias/{condition}-v-{control}_tss-seq-spikenorm-intragenic-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down"], reference=["ATG", "RELATIVE", "STOP"]),
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/antisense/position_bias/{condition}-v-{control}_tss-seq-libsizenorm-antisense-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down"], reference=["TSS", "RELATIVE", "CPS"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/antisense/position_bias/{condition}-v-{control}_tss-seq-spikenorm-antisense-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down"], reference=["TSS", "RELATIVE", "CPS"]),
        ##distances of differentially expressed peaks
        #expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{{ttype}}/{condition}-v-{control}-relative-distances-libsizenorm-{{direction}}-{{ttype}}.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], ttype=["intragenic", "antisense"]),
        #expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{{ttype}}/{condition}-v-{control}-relative-distances-spikenorm-{{direction}}-{{ttype}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], ttype=["intragenic", "antisense"]),
        #datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/tss-seq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup-sense.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, status=["all","passing"]) if config["plot_figures"] else [],
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/tss-seq_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup-sense.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, status=["all","passing"]) if config["plot_figures"] else [],
        #correlations of transcript classes with genic TSSs
        expand("diff_exp/{condition}-v-{control}/libsizenorm/class_v_genic/{condition}-v-{control}_tss-seq-libsizenorm-class-v-genic.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/spikenorm/class_v_genic/{condition}-v-{control}_tss-seq-spikenorm-class-v-genic.tsv", zip, condition=conditiongroups_si, control=controlgroups_si),
        #enrichment of known motifs
        expand(expand("motifs/{condition}-v-{control}/libsizenorm/{{category}}/{{negative}}/{condition}-v-{control}_tss-seq-libsizenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES) if config["motifs"]["run_motif_analyses"] else [],
        expand(expand("motifs/{condition}-v-{control}/spikenorm/{{category}}/{{negative}}/{condition}-v-{control}_tss-seq-spikenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES) if config["motifs"]["run_motif_analyses"] else [],
        ## MEME-ChIP
        #expand(expand("motifs/meme/{condition}-v-{control}/libsizenorm/{{region}}/{condition}-v-{control}-results-libsizenorm-{{direction}}-{{category}}-{{region}}-meme_chip/meme-chip.html", zip, condition=conditiongroups, control=controlgroups), region=["upstream","peak"], direction=["up", "down"], category=CATEGORIES) if config["motifs"]["meme-chip"]["run-meme-chip"] else [],
        #expand(expand("motifs/meme/{condition}-v-{control}/spikenorm/{{region}}/{condition}-v-{control}-results-spikenorm-{{direction}}-{{category}}-{{region}}-meme_chip/meme-chip.html", zip, condition=conditiongroups_si, control=controlgroups_si), region=["upstream","peak"], direction=["up", "down"], category=CATEGORIES) if config["motifs"]["meme-chip"]["run-meme-chip"] else [],
        #gene ontology
        expand(expand("gene_ontology/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_tss-seq-libsizenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down", "unchanged"], category=["genic", "intragenic", "antisense", "convergent", "divergent"]) if config["run_gene_ontology"] else [],
        expand(expand("gene_ontology/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_tss-seq-spikenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down", "unchanged"], category=["genic", "intragenic", "antisense", "convergent", "divergent"]) if config["run_gene_ontology"] else [],
        #sequence logos
        expand("seq_logos/{group}/{group}-seqlogos.svg", group=set([PASSING[x]['group'] for x in PASSING])),

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

rule intragenic_position_bias:
    input:
        diffexp_results = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}.tsv",
        chrom_sizes = config["genome"]["chrsizes"],
        genic_regions = build_annotations(os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
    output:
        atg = "diff_exp/{condition}-v-{control}/{norm}/intragenic/position_bias/{condition}-v-{control}_tss-seq-{norm}-intragenic-{direction}-ATG.tsv",
        rel = "diff_exp/{condition}-v-{control}/{norm}/intragenic/position_bias/{condition}-v-{control}_tss-seq-{norm}-intragenic-{direction}-RELATIVE.tsv",
        stop = "diff_exp/{condition}-v-{control}/{norm}/intragenic/position_bias/{condition}-v-{control}_tss-seq-{norm}-intragenic-{direction}-STOP.tsv",
    log: "logs/intragenic_position_bias/intragenic_position_bias-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/intragenic_tss_position_bias.py -d {input.diffexp_results} -c {input.chrom_sizes} -g {input.genic_regions} -a {output.atg} -r {output.rel} -s {output.stop}) &> {log}
        """

rule antisense_position_bias:
    input:
        diffexp_results = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-antisense-{direction}.tsv",
        chrom_sizes = config["genome"]["chrsizes"],
    output:
        tss = "diff_exp/{condition}-v-{control}/{norm}/antisense/position_bias/{condition}-v-{control}_tss-seq-{norm}-antisense-{direction}-TSS.tsv",
        rel = "diff_exp/{condition}-v-{control}/{norm}/antisense/position_bias/{condition}-v-{control}_tss-seq-{norm}-antisense-{direction}-RELATIVE.tsv",
        cps = "diff_exp/{condition}-v-{control}/{norm}/antisense/position_bias/{condition}-v-{control}_tss-seq-{norm}-antisense-{direction}-CPS.tsv",
    log: "logs/antisense_position_bias/antisense_position_bias-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/antisense_tss_position_bias.py -d {input.diffexp_results} -c {input.chrom_sizes} -a {output.tss} -r {output.rel} -s {output.cps}) &> {log}
        """

rule genic_v_class:
    input:
        genic = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-all.tsv",
        intragenic = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-all.tsv",
        antisense = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-antisense-all.tsv",
        convergent = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-all.tsv",
        divergent = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-all.tsv",
    output:
        tsv = "diff_exp/{condition}-v-{control}/{norm}/class_v_genic/{condition}-v-{control}_tss-seq-{norm}-class-v-genic.tsv",
        lfc_v_lfc = "diff_exp/{condition}-v-{control}/{norm}/class_v_genic/{condition}-v-{control}_tss-seq-{norm}-class-v-genic-lfc-v-lfc.svg",
        lfc_v_expr = "diff_exp/{condition}-v-{control}/{norm}/class_v_genic/{condition}-v-{control}_tss-seq-{norm}-class-v-genic-lfc-v-expr.svg",
        expr_v_expr = "diff_exp/{condition}-v-{control}/{norm}/class_v_genic/{condition}-v-{control}_tss-seq-{norm}-class-v-genic-expr-v-expr.svg",
    script: "scripts/tss_class_v_genic.R"

