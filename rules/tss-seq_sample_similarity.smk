#!/usr/bin/env python

rule union_bedgraph:
    input:
        expand("coverage/{{norm}}/{sample}_tss-seq-{{norm}}-SENSE.bedgraph", sample=SAMPLES)
    output:
        "qual_ctrl/scatter_plots/tss-seq_union-bedgraph-{norm}-allsamples.tsv.gz"
    params:
        names = list(SAMPLES.keys())
    log: "logs/union_bedgraph-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz > {output}) &> {log}
        """

#TODO: plot correlations over multiple binsizes, and over annotations
# in the current implementation the bins are not actually guaranteed to be a single base,
# although in practice they are
rule plot_scatter_plots:
    input:
        "qual_ctrl/scatter_plots/tss-seq_union-bedgraph-{norm}-allsamples.tsv.gz"
    output:
        "qual_ctrl/scatter_plots/{condition}-v-{control}/{status}/{condition}-v-{control}_tss-seq-{norm}-scatterplots-{status}.svg"
    params:
        pcount = 0.1,
        samplelist = get_condition_control_samples
    script:
        "../scripts/plot_scatter_plots.R"

