#!/usr/bin/env python

import os
import re
import subprocess
import itertools
from math import log2, log10

configfile: "config.yaml"

subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

configfile: build_annotations("config.yaml")

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}

#groups which have at least two passing samples, so that they are valid for peakcalling and diff exp
validgroups = set(z for z in [PASSING[x]['group'] for x in PASSING] if [PASSING[x]['group'] for x in PASSING].count(z)>=2)
validgroups_si = set(z for z in [PASSING[x]['group'] for x in PASSING if PASSING[x]['spikein']] if [PASSING[x]['group'] for x in PASSING].count(z)>=2)

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))

comparisons_si = config["comparisons"]["spikenorm"]
if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in config["comparisons"]["spikenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups_si]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in config["comparisons"]["spikenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups_si]))

CATEGORIES = ["genic", "intragenic", "antisense", "convergent", "divergent", "intergenic"]

FIGURES = config["figures"]

#get all motif names from motif databases, cleaning nasty characters in some motif names
MOTIFS = set(subprocess.run(args="meme2meme " + " ".join(config["motifs"]["databases"]) + " | grep -e '^MOTIF' | cut -d ' ' -f2 | sed 's/\//_/g; s/&/_/g; s/{/[/g; s/}/]/g' ", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()) if config["motifs"]["run_motif_analyses"] else []

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys())),
    group = "|".join(set(re.escape(v["group"]) for k,v in SAMPLES.items())),
    control = "|".join(set(re.escape(x) for x in controlgroups + ([] if not comparisons_si else controlgroups_si) + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + ([] if not comparisons_si else conditiongroups_si) + ["all"])),
    species = "experimental|spikein",
    read_status = "raw|cleaned|aligned_noPCRdup|unaligned",
    category = "|".join(CATEGORIES + ["all"]),
    figure = "|".join(re.escape(x) for x in list(FIGURES.keys())),
    annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES]))),
    motif = "|".join(re.escape(x) for x in MOTIFS),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|sicounts|libsizenorm|spikenorm",
    strand = "SENSE|ANTISENSE|plus|minus",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    negative = "unchanged|random"

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
    genic_v_class,

def statuscheck(dict1, dict2):
    return(["passing"] if dict1 == dict2 else ["all", "passing"])

def conditioncheck(conditionlist):
    return(conditionlist if len(conditionlist)==1 else conditionlist + ["all"])

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #FastQC
        'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bw", sample=SAMPLES, norm=["counts", "libsizenorm"], strand=["SENSE","ANTISENSE","plus","minus"]),
        expand("coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bw", sample=SISAMPLES, norm=["sicounts", "spikenorm"], strand=["SENSE","ANTISENSE","plus","minus"]),
        #quality controls
        "qual_ctrl/read_processing/tss-seq_read-processing-loss.svg",
        expand("qual_ctrl/spikein/tss-seq_spikein-plots-{status}.svg", status=statuscheck(SISAMPLES, SIPASSING)) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_tss-seq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), status=statuscheck(SAMPLES, PASSING), windowsize = config["scatterplot_binsizes"]),
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_tss-seq-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), status=statuscheck(SISAMPLES, SIPASSING), windowsize = config["scatterplot_binsizes"]) if SISAMPLES and comparisons_si else [],
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
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_spikein-peaks.bed", zip, condition=conditiongroups_si, control=controlgroups_si) if comparisons_si else [],
        #differential expression of peaks
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-experimental-peak-counts.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/{condition}-v-{control}_tss-seq-spikein-peak-counts.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if comparisons_si else [],
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-results-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups),direction=["all", "up", "unchanged", "down"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{condition}-v-{control}_tss-seq-spikenorm-diffexp-results-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si),direction=["all", "up", "unchanged", "down"]) if comparisons_si else [],
        #categorize differentially expressed peaks
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-results-{{category}}-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups), direction = ["all","up","unchanged","down"], category=CATEGORIES),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_tss-seq-spikenorm-diffexp-results-{{category}}-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si), direction = ["all","up","unchanged","down"], category=CATEGORIES) if comparisons_si else [],
        #differential expression summary
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups, control=controlgroups), plot = ["mosaic", "maplot", "volcano"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{condition}-v-{control}_tss-seq-spikenorm-diffexp-{{plot}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), plot = ["mosaic", "maplot", "volcano"]) if comparisons_si else [],
        #random distributions for differentially expressed peaks
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/intragenic/position_bias/{condition}-v-{control}_tss-seq-libsizenorm-intragenic-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down"], reference=["ATG", "RELATIVE", "STOP"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/intragenic/position_bias/{condition}-v-{control}_tss-seq-spikenorm-intragenic-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down"], reference=["ATG", "RELATIVE", "STOP"]),
        expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/antisense/position_bias/{condition}-v-{control}_tss-seq-libsizenorm-antisense-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down"], reference=["TSS", "RELATIVE", "CPS"]),
        expand(expand("diff_exp/{condition}-v-{control}/spikenorm/antisense/position_bias/{condition}-v-{control}_tss-seq-spikenorm-antisense-{{direction}}-{{reference}}.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down"], reference=["TSS", "RELATIVE", "CPS"]),
        #distances of differentially expressed peaks
        #expand(expand("diff_exp/{condition}-v-{control}/libsizenorm/{{ttype}}/{condition}-v-{control}-relative-distances-libsizenorm-{{direction}}-{{ttype}}.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], ttype=["intragenic", "antisense"]),
        #expand(expand("diff_exp/{condition}-v-{control}/spikenorm/{{ttype}}/{condition}-v-{control}-relative-distances-spikenorm-{{direction}}-{{ttype}}.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], ttype=["intragenic", "antisense"]),
        #datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/tss-seq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup-sense.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), figure=FIGURES, status=statuscheck(SAMPLES, PASSING)) if config["plot_figures"] else [],
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/tss-seq_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup-sense.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), figure=FIGURES, status=statuscheck(SISAMPLES, SIPASSING)) if config["plot_figures"] and SISAMPLES and comparisons_si else [],
        #correlations of transcript classes with genic TSSs
        expand("diff_exp/{condition}-v-{control}/libsizenorm/class_v_genic/{condition}-v-{control}_tss-seq-libsizenorm-class-v-genic.tsv", zip, condition=conditiongroups, control=controlgroups),
        expand("diff_exp/{condition}-v-{control}/spikenorm/class_v_genic/{condition}-v-{control}_tss-seq-spikenorm-class-v-genic.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if comparisons_si else[],
        #enrichment of known motifs
        expand(expand("motifs/{condition}-v-{control}/libsizenorm/{{category}}/{{negative}}/{condition}-v-{control}_tss-seq-libsizenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES) if config["motifs"]["run_motif_analyses"] else [],
        expand(expand("motifs/{condition}-v-{control}/spikenorm/{{category}}/{{negative}}/{condition}-v-{control}_tss-seq-spikenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES) if config["motifs"]["run_motif_analyses"] and comparisons_si else [],
        #gene ontology
        expand(expand("gene_ontology/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_tss-seq-libsizenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down", "unchanged"], category=["genic", "intragenic", "antisense", "convergent", "divergent"]) if config["run_gene_ontology"] else [],
        expand(expand("gene_ontology/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_tss-seq-spikenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down", "unchanged"], category=["genic", "intragenic", "antisense", "convergent", "divergent"]) if config["run_gene_ontology"] and comparisons_si else [],
        #sequence logos
        expand("seq_logos/{group}/{group}-seqlogos.svg", group=set([PASSING[x]['group'] for x in PASSING])),
        expand("seq_logos/{group}/{group}-{category}-seqlogo.meme", group=set([PASSING[x]['group'] for x in PASSING]), category=CATEGORIES+["all"]),
        # MEME-ChIP
        expand(expand("motifs/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_tss-seq-libsizenorm-diffexp-results-{{category}}-{{direction}}-meme_chip/summary.tsv", zip, condition=conditiongroups, control=controlgroups), category=CATEGORIES, direction=["up", "down", "unchanged"]) if config["motifs"]["meme-chip"]["run-meme-chip"] else [],
        expand(expand("motifs/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_tss-seq-spikenorm-diffexp-results-{{category}}-{{direction}}-meme_chip/summary.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), category=CATEGORIES, direction=["up", "down", "unchanged"]) if config["motifs"]["meme-chip"]["run-meme-chip"] and comparisons_si else [],

rule intragenic_position_bias:
    input:
        diffexp_results = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}.tsv",
        fasta = config["genome"]["fasta"],
        genic_regions = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
    output:
        atg = "diff_exp/{condition}-v-{control}/{norm}/intragenic/position_bias/{condition}-v-{control}_tss-seq-{norm}-intragenic-{direction}-ATG.tsv",
        rel = "diff_exp/{condition}-v-{control}/{norm}/intragenic/position_bias/{condition}-v-{control}_tss-seq-{norm}-intragenic-{direction}-RELATIVE.tsv",
        stop = "diff_exp/{condition}-v-{control}/{norm}/intragenic/position_bias/{condition}-v-{control}_tss-seq-{norm}-intragenic-{direction}-STOP.tsv",
    conda: "envs/position_bias.yaml"
    log: "logs/intragenic_position_bias/intragenic_position_bias-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/intragenic_tss_position_bias.py -d {input.diffexp_results} -c $(faidx {input.fasta} -i chromsizes) -g {input.genic_regions} -a {output.atg} -r {output.rel} -s {output.stop}) &> {log}
        """

rule antisense_position_bias:
    input:
        diffexp_results = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-antisense-{direction}.tsv",
        fasta = config["genome"]["fasta"],
    output:
        tss = "diff_exp/{condition}-v-{control}/{norm}/antisense/position_bias/{condition}-v-{control}_tss-seq-{norm}-antisense-{direction}-TSS.tsv",
        rel = "diff_exp/{condition}-v-{control}/{norm}/antisense/position_bias/{condition}-v-{control}_tss-seq-{norm}-antisense-{direction}-RELATIVE.tsv",
        cps = "diff_exp/{condition}-v-{control}/{norm}/antisense/position_bias/{condition}-v-{control}_tss-seq-{norm}-antisense-{direction}-CPS.tsv",
    conda: "envs/position_bias.yaml"
    log: "logs/antisense_position_bias/antisense_position_bias-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (python scripts/antisense_tss_position_bias.py -d {input.diffexp_results} -c $(faidx {input.fasta} -i chromsizes) -a {output.tss} -r {output.rel} -s {output.cps}) &> {log}
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
    conda: "envs/tidyverse.yaml"
    script: "scripts/tss_class_v_genic.R"

