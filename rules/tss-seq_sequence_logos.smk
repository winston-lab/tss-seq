#!/usr/bin/env python

localrules: plot_seqlogos

rule get_seqlogo_data:
    input:
        bam = lambda wc: expand("alignment/{sample}_tss-seq-noPCRduplicates.bam", sample=getsamples(wc.group, wc.group)),
        bed = lambda wc: "peakcalling/{group}/{category}/{group}-experimental-idrpeaks-{category}.narrowpeak".format(**wc) if wc.category != "all" else [],
        chrsizes = config["genome"]["chrsizes"],
        fasta = config["genome"]["fasta"]
    params:
        prefix = config["combinedgenome"]["experimental_prefix"],
        slop = int(config["consensus"]["window"]),
        windowsize = int(2*config["consensus"]["window"]+1),
        intersect_cmd = lambda wc: "bedtools intersect -wa -s -a stdin -b peakcalling/{group}/{category}/{group}-experimental-idrpeaks-{category}.narrowpeak | ".format(**wc) if wc.category != "all" else ""
    output:
        seqlogo_data = "seq_logos/{group}/{group}-{category}-seqlogo.tsv",
    log: "logs/get_seqlogo_data/get_seqlogo_data-{group}-{category}.log"
    run:
        fasta_stats = subprocess.run(args="fasta-get-markov " + input.fasta + " | tail -n +4", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()
        composition = "{" + ",".join(["'" + fasta_stats[x] + "':" + str(1e2*float(fasta_stats[x+1])) for x in range(0,8,2)]) + "}"
        shell("""(samtools merge - {input.bam} | bedtools bamtobed -i stdin | sed -n 's/{params.prefix}//p' | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6}}' | {params.intersect_cmd} bedtools slop -b {params.slop} -s -i stdin -g {input.chrsizes} | awk '$3-$2=={params.windowsize}' | bedtools getfasta -s -fi {input.fasta} -bed stdin | weblogo --format logodata --sequence-type dna --composition "{composition}" > {output.seqlogo_data}) &> {log}""")

rule plot_seqlogos:
    input:
        seqlogo_data = expand("seq_logos/{{group}}/{{group}}-{category}-seqlogo.tsv", category = ["all"] + CATEGORIES),
        fasta = config["genome"]["fasta"]
    params:
        tss_classes = ["all"] + CATEGORIES,
        slop = int(config["consensus"]["window"]),
    output:
        "seq_logos/{group}/{group}-seqlogos.svg",
    log: "logs/plot_seqlogos/plot_seqlogos-{group}.log"
    run:
        gc_pct = float(subprocess.run(args="fasta-get-markov " + input.fasta + " | tail -n +4", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()[3])*2
        shell(""" (Rscript scripts/plot_seqlogo.R --input {input.seqlogo_data} --tss_classes {params.tss_classes} -b {params.slop} --gc_pct {gc_pct} -l {wildcards.group} -o {output}) &> {log}""")

