#!/usr/bin/env python

localrules: plot_seqlogos, seqlogo_to_meme

rule get_seqlogo_data:
    input:
        bam = lambda wc: expand("alignment/{sample}_tss-seq-noPCRduplicates-experimental.bam", sample=get_samples("passing", "libsizenorm", wc.group)),
        bed = lambda wc: "peakcalling/{group}/{category}/{group}-experimental-idrpeaks-{category}.narrowpeak".format(**wc) if wc.category != "all" else [],
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        seqlogo_data = "seq_logos/{group}/{group}-{category}-seqlogo.tsv",
    params:
        slop = int(config["consensus"]["window"]),
        windowsize = int(2*config["consensus"]["window"]+1),
        intersect_cmd = lambda wc: "bedtools intersect -wa -s -a stdin -b peakcalling/{group}/{category}/{group}-experimental-idrpeaks-{category}.narrowpeak | ".format(**wc) if wc.category != "all" else ""
    log:
        "logs/get_seqlogo_data/get_seqlogo_data-{group}-{category}.log"
    run:
        fasta_stats = subprocess.run(args="fasta-get-markov " + input.fasta + " | tail -n +4", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()
        composition = "{" + ",".join(["'" + fasta_stats[x] + "':" + str(1e2*float(fasta_stats[x+1])) for x in range(0,8,2)]) + "}"
        shell("""(samtools merge - {input.bam} | \
                  bedtools bamtobed -i stdin | \
                  awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6}}' | \
                  {params.intersect_cmd} bedtools slop -b {params.slop} -s -i stdin -g <(faidx {input.fasta} -i chromsizes) | \
                  awk '$3-$2=={params.windowsize}' | \
                  bedtools getfasta -s -fi {input.fasta} -bed stdin | \
                  weblogo --format logodata --sequence-type dna --composition "{composition}" > {output.seqlogo_data}) &> {log}""")

rule plot_seqlogos:
    input:
        seqlogo_data = expand("seq_logos/{{group}}/{{group}}-{category}-seqlogo.tsv", category = ["all"] + CATEGORIES),
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    params:
        tss_classes = ["all"] + CATEGORIES,
        slop = int(config["consensus"]["window"]),
    output:
        "seq_logos/{group}/{group}-seqlogos.svg",
    log:
        "logs/plot_seqlogos/plot_seqlogos-{group}.log"
    run:
        gc_pct = float(subprocess.run(args="fasta-get-markov " + input.fasta + " | tail -n +4", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()[3])*2
        shell(""" (Rscript scripts/plot_seqlogo.R --input {input.seqlogo_data} --tss_classes {params.tss_classes} -b {params.slop} --gc_pct {gc_pct} -l {wildcards.group} -o {output}) &> {log}""")

rule seqlogo_to_meme:
    input:
        data = "seq_logos/{group}/{group}-{category}-seqlogo.tsv",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        "seq_logos/{group}/{group}-{category}-seqlogo.meme",
    log:
        "logs/seqlogo_to_meme/seqlogo_to_meme-{group}-{category}.log"
    shell: """
        (grep -v -e "#" {input.data} | cut -f2-5 | matrix2meme -numseqs $(grep -v -e "#" {input.data} | head -n 1 | cut -f2-5 | awk '{{print $1+$2+$3+$4}}') -logodds -bg <(fasta-get-markov {input.fasta}) > {output}) &> {log}
        """

