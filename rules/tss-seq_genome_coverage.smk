#!/usr/bin/env python

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
    log: "logs/genome_coverage/genome_coverage_{sample}-{counttype}-{strand}.log"
    shell: """
        (bedtools genomecov -bga -5 -strand {params.strand} -ibam {input} | LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

#NOTE: although we could do this by looking up library size values
# from the library sizes file, this way we do not have to wait for
# all other samples finish aligning in order to normalize
rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_tss-seq-counts-{strand}.bedgraph",
        bam = lambda wc: f"alignment/{wc.sample}_tss-seq-noPCRduplicates-" + ("experimental" if wc.norm=="spikenorm" else "spikein") + ".bam",
    output:
        normalized = "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bedgraph",
    params:
        scale_factor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1,
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus"
    log: "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{strand}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam} | paste -d "" - <(echo "/({params.scale_factor}*1000000)") | bc -l) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = "coverage/{norm}/{sample}_tss-seq-{norm}-plus.bedgraph",
        minus = "coverage/{norm}/{sample}_tss-seq-{norm}-minus.bedgraph"
    output:
        sense = "coverage/{norm}/{sample}_tss-seq-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}_tss-seq-{norm}-ANTISENSE.bedgraph"
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus} > {output.antisense}) &>> {log}
        """

def select_chromsizes(wc):
    if wc.strand in ["plus", "minus"]:
        if wc.norm=="sicounts":
            return config["genome"]["sichrsizes"]
        return config["genome"]["chrsizes"]
    if wc.norm=="sicounts":
        return os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv"
    return os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bedgraph",
        chrsizes = select_chromsizes
    output:
        "coverage/{norm}/{sample}_tss-seq-{norm}-{strand}.bw",
    log : "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

