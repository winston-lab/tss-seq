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
    script: "../scripts/processing_summary.R"

rule build_spikein_counts_table:
    input:
        bams = expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=sisamples),
        bais = expand("alignment/{sample}_tss-seq-noPCRduplicates.bam.bai", sample=sisamples),
        chrsizes = config["combinedgenome"]["chrsizes"]
    output:
        "qual_ctrl/spikein/tss-seq_spikein-counts.tsv"
    params:
        groups = [v["group"] for k,v in sisamples.items()],
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
        spikein_prefix = config["combinedgenome"]["spikein_prefix"],
    log: "logs/build_spikein_counts_table.log"
    run:
        shell("""(echo -e "sample\tgroup\ttotal_counts\texperimental_counts\tspikein_counts" > {output}) &> {log} """)
        for sample, group, bam in zip(sisamples.keys(), params.groups, input.bams):
            shell("""(paste <(echo -e "{sample}\t{group}") <(samtools view -c {bam}) <(grep -oh "\w*{params.exp_prefix}\w*" {input.chrsizes} | xargs samtools view -c {bam}) <(grep -oh "\w*{params.spikein_prefix}\w*" {input.chrsizes} | xargs samtools view -c {bam}) >> {output}) &>> {log}""")

rule plot_spikein_pct:
    input:
        "qual_ctrl/spikein/tss-seq_spikein-counts.tsv"
    output:
        plot = "qual_ctrl/spikein/tss-seq_spikein-plots-{status}.svg",
        stats = "qual_ctrl/spikein/tss-seq_spikein-stats-{status}.tsv"
    params:
        samplelist = lambda wc : list(sisamples.keys()) if wc.status=="all" else list(sipassing.keys()),
        conditions = conditiongroups_si,
        controls = controlgroups_si,
    script: "../scripts/plot_si_pct.R"

