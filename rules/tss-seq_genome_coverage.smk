#!/usr/bin/env python

localrules:
    normalize_genome_coverage,
    make_stranded_bedgraph,
    bedgraph_to_bigwig

rule genome_coverage:
    input:
        lambda wc: f"alignment/{wc.sample}_tss-seq-noPCRduplicates-" + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam"
    output:
        "coverage/{counttype}/{sample}_tss-seq-{counttype}-{strand}.bedgraph",
    params:
        strand = lambda wc: {"plus": "+", "minus": "-"}.get(wc.strand)
    wildcard_constraints:
        counttype="counts|sicounts",
        strand="plus|minus"
    log:
        "logs/genome_coverage/genome_coverage_{sample}-{counttype}-{strand}.log"
    shell: """
        (bedtools genomecov -bga -5 -strand {params.strand} -ibam {input} | \
         LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

#NOTE: although we could do this by looking up library size values
# from the library sizes file, this way we do not have to wait for
# all other samples finish aligning in order to normalize
rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_tss-seq-counts-{strand}.bedgraph",
        bam = lambda wc: f"alignment/{wc.sample}_tss-seq-noPCRduplicates-" + ("spikein" if wc.norm=="spikenorm" else "experimental") + ".bam",
    output:
        normalized = "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bedgraph",
    params:
        scale_factor = lambda wc: config["spike_in"]["proportion"] if wc.norm=="spikenorm" else 1,
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus"
    log:
        "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{strand}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam} | paste -d "" - <(echo "/({params.scale_factor}*1000000)") | bc -l) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = lambda wc: "coverage/{norm}/{sample}_tss-seq-{norm}-".format(**wc) + ("plus" if wc.strand=="SENSE" else "minus") + ".bedgraph",
        minus = lambda wc: "coverage/{norm}/{sample}_tss-seq-{norm}-".format(**wc) + ("minus" if wc.strand=="SENSE" else "plus") + ".bedgraph",
    output:
        "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bedgraph",
    wildcard_constraints:
        strand="SENSE|ANTISENSE"
    log:
        "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}-{strand}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bedgraph",
        fasta = lambda wc: os.path.abspath(config["spike_in"]["fasta"]) if wc.norm=="sicounts" else os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bw",
    params:
        stranded = lambda wc: [] if wc.strand in ["plus", "minus"] else """| awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2; print $1"-minus", $2}}' | LC_COLLATE=C sort -k1,1"""
    log:
        "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes {params.stranded}) {output}) &> {log}
        """

