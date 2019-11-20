#!/usr/bin/env python

localrules:
    make_stranded_annotations,
    compute_matrix,
    cat_matrices

rule make_stranded_annotations:
    input:
        lambda wc : FIGURES[wc.figure]["annotations"][wc.annotation]["path"]
    output:
        "datavis/{figure}/{annotation}.bed"
    log:
        "logs/make_stranded_annotations/make_stranded_annotations-{figure}-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

rule compute_matrix:
    input:
        annotation = "datavis/{figure}/{annotation}.bed",
        bw = "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bw"
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
    threads:
        config["threads"]
    log:
        "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}-{norm}-{strand}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} -g {params.group} -s {wildcards.sample} -a {params.anno_label} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}-melted.tsv.gz", annotation=list(FIGURES[wc.figure]["annotations"].keys()), sample=SAMPLES, figure=wc.figure, norm=wc.norm, strand=wc.strand)
    output:
        "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-tss-seq-{norm}-{strand}.tsv.gz"
    log:
        "logs/cat_matrices/cat_matrices-{figure}_{norm}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        matrices = expand("datavis/{{figure}}/{{norm}}/{{figure}}-allsamples-allannotations-tss-seq-{{norm}}-{strand}.tsv.gz", strand=["SENSE", "ANTISENSE"]),
        annotations = lambda wc: [v["path"] for k,v in FIGURES[wc.figure]["annotations"].items()]
    output:
        heatmap_sample_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample-bothstrands.svg",
        heatmap_sample_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample-sense.svg",
        heatmap_sample_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample-antisense.svg",
        heatmap_group_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup-bothstrands.svg",
        heatmap_group_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup-sense.svg",
        heatmap_group_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup-antisense.svg",
        metagene_sample_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-bothstrands.svg",
        metagene_sample_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-sense.svg",
        metagene_sample_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-antisense.svg",
        metagene_sample_overlay_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-bothstrands.svg",
        metagene_sample_overlay_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-sense.svg",
        metagene_sample_overlay_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-antisense.svg",
        metagene_group_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-bothstrands.svg",
        metagene_group_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-sense.svg",
        metagene_group_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-antisense.svg",
        metagene_sampleclust_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-bothstrands.svg",
        metagene_sampleclust_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-sense.svg",
        metagene_sampleclust_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-antisense.svg",
        metagene_groupclust_both = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-bothstrands.svg",
        metagene_groupclust_sense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-sense.svg",
        metagene_groupclust_antisense = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/tss-seq_{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-antisense.svg",
    params:
        # abusing snakemake a bit here...using params as output paths in order to use lambda functions
        annotations_out = lambda wc: ["datavis/{figure}/{norm}/{condition}-v-{control}/{status}/".format(**wc) + annotation + "_cluster-" + str(cluster) + ".bed" for annotation in FIGURES[wc.figure]["annotations"] for cluster in range(1, FIGURES[wc.figure]["annotations"][annotation]["n_clusters"]+1)],
        clusters_out = lambda wc: ["datavis/{figure}/{norm}/{condition}-v-{control}/{status}/".format(**wc) + annotation + ".pdf" for annotation in FIGURES[wc.figure]["annotations"]],
        samplelist = lambda wc: get_samples(wc.status, wc.norm, [wc.condition, wc.control]),
        plottype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        pct_cutoff = lambda wc: FIGURES[wc.figure]["parameters"]["pct_cutoff"],
        log_transform = lambda wc: str(FIGURES[wc.figure]["parameters"]["log_transform"]).upper(),
        pcount = lambda wc: 0 if not FIGURES[wc.figure]["parameters"]["log_transform"] else FIGURES[wc.figure]["parameters"]["pseudocount"],
        spread_type = lambda wc: FIGURES[wc.figure]["parameters"]["spread_type"],
        trim_pct = lambda wc: FIGURES[wc.figure]["parameters"]["trim_pct"],
        refpointlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
        endlabel = lambda wc:  "HAIL SATAN" if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["three_prime_label"],
        cmap = lambda wc: FIGURES[wc.figure]["parameters"]["heatmap_colormap"],
        sortmethod = lambda wc: FIGURES[wc.figure]["parameters"]["arrange"],
        cluster_scale = lambda wc: "FALSE" if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else str(FIGURES[wc.figure]["parameters"]["cluster_scale"]).upper(),
        cluster_samples = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else cluster_samples(wc.status, wc.norm, list(FIGURES[wc.figure]["parameters"]["cluster_conditions"].keys()), list(FIGURES[wc.figure]["parameters"]["cluster_conditions"].values())),
        cluster_five = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_five"],
        cluster_three = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_three"],
        k = lambda wc: [v["n_clusters"] for k,v in FIGURES[wc.figure]["annotations"].items()],
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_tss_figures.R"

