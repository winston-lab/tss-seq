#!/usr/bin/env python

rule gene_ontology:
    input:
        universe = lambda wc: config["genome"]["orf_annotation"] if wc.category=="intragenic" else config["genome"]["transcript_annotation"],
        diffexp_path = "diff_exp/peaks/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-{category}-{direction}.tsv",
        go_anno_path = config["gene_ontology_mapping_file"]
    output:
        results = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-gene-ontology-results.tsv",
        enriched_combined = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-gene-ontology-enriched-all.svg",
        depleted_combined = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-gene-ontology-depleted-all.svg",
        enriched_facet = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-gene-ontology-enriched-facetted.svg",
        depleted_facet = "gene_ontology/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-{category}-{direction}-gene-ontology-depleted-facetted.svg",
    params:
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated" if wc.direction=="down" else wc.direction
    conda:
        "../envs/gene_ontology.yaml"
    script:
        "../scripts/gene_ontology.R"

