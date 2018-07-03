#!/usr/bin/env python

localrules: aggregate_read_numbers,
    build_spikein_counts_table,
    plot_spikein_pct

rule aggregate_read_numbers:
    input:
        adapter = expand("logs/clean_reads/remove_adapter-{sample}.log", sample=SAMPLES),
        qual_trim = expand("logs/clean_reads/remove_3p_bc_and_trim-{sample}.log", sample=SAMPLES),
        align = expand("alignment/{sample}/align_summary.txt", sample=SAMPLES),
        nodups = expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=SAMPLES)
    output:
        "qual_ctrl/read_processing/tss-seq_read-processing-summary.tsv"
    log: "logs/aggregate_read_numbers.log"
    run:
        shell("""(echo -e "sample\traw\tadapter_removed\tquality_trimmed\tmapped\tunique_map\tnoPCRdup" > {output}) &> {log}""")
        for sample, adapter, qual_trim, align, nodups in zip(SAMPLES.keys(), input.adapter, input.qual_trim, input.align, input.nodups):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &>> {log}""")
            shell("""(grep -e "Reads written" {qual_trim} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"}}{{print $1}}' >> {output}) &>> {log}""")
            shell("""(awk 'BEGIN{{ORS="\t"}} NR==3 || NR==4{{print $3}}' {align} >> {output}) &> {log}""")
            shell("""(samtools view -c {nodups} | awk '{{print $1}}' >> {output}) &>> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$6=$5-$6; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &>> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing/tss-seq_read-processing-summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing/tss-seq_read-processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing/tss-seq_read-processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing/tss-seq_read-processing-loss.svg",
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/processing_summary.R"

rule build_spikein_counts_table:
    input:
        total_bam = expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=SISAMPLES),
        exp_bam = expand("alignment/{sample}_tss-seq-noPCRduplicates-experimental.bam", sample=SISAMPLES),
        si_bam = expand("alignment/{sample}_tss-seq-noPCRduplicates-spikein.bam", sample=SISAMPLES),
    output:
        "qual_ctrl/spikein/tss-seq_spikein-counts.tsv"
    params:
        groups = [v["group"] for k,v in SISAMPLES.items()],
    log: "logs/build_spikein_counts_table.log"
    run:
        shell("""(echo -e "sample\tgroup\ttotal_counts\texperimental_counts\tspikein_counts" > {output}) &> {log} """)
        for sample, group, total, exp, si in zip(SISAMPLES.keys(), params.groups, input.total_bam, input.exp_bam, input.si_bam):
            shell("""(paste <(echo -e "{sample}\t{group}") <(samtools view -c {total}) <(samtools view -c {exp}) <(samtools view -c {si}) >> {output}) &>> {log}""")

rule plot_spikein_pct:
    input:
        "qual_ctrl/spikein/tss-seq_spikein-counts.tsv"
    output:
        plot = "qual_ctrl/spikein/tss-seq_spikein-plots-{status}.svg",
        stats = "qual_ctrl/spikein/tss-seq_spikein-stats-{status}.tsv"
    params:
        samplelist = lambda wc : list(SISAMPLES.keys()) if wc.status=="all" else list(SIPASSING.keys()),
        conditions = [] if not SISAMPLES else conditiongroups_si,
        controls = [] if not SISAMPLES else controlgroups_si,
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/plot_si_pct.R"

