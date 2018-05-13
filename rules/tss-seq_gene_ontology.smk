#!/usr/bin/env python

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
