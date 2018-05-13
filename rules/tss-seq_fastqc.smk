#!/usr/bin/env python

#fastqc for raw or cleaned reads
rule fastqc_prealignment:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"] if wc.read_status=="raw" else "fastq/cleaned/" + wc.sample + "_tss-seq-clean.fastq.gz"
    output:
        "qual_ctrl/fastqc/{read_status}/{sample}_fastqc-data-{read_status}.txt"
    params:
        fname = lambda wc: re.split('.fq|.fastq', os.path.split(SAMPLES[wc.sample]["fastq"])[1])[0] if wc.read_status=="raw" else wc.sample + "_tss-seq-clean",
        adapter = config["cutadapt"]["adapter"]
    wildcard_constraints:
        read_status="raw|cleaned"
    threads : config["threads"]
    log: "logs/fastqc/fastqc_{read_status}-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.read_status}) &> {log}
        (fastqc --adapters <(echo -e "adapter\t{params.adapter}") --nogroup --noextract -t {threads} -o qual_ctrl/fastqc/{wildcards.read_status} {input}) &>> {log}
        (unzip -p qual_ctrl/fastqc/{wildcards.read_status}/{params.fname}_fastqc.zip {params.fname}_fastqc/fastqc_data.txt > {output}) &>> {log}
        """

#fastqc for no PCR duplicates and unalignd reads
rule fastqc_postalignment:
    input:
        lambda wc: "alignment/" + wc.sample + "_tss-seq-noPCRduplicates.bam" if wc.read_status=="aligned_noPCRdup" else "alignment/" + wc.sample + "/unmapped.bam",
    output:
        "qual_ctrl/fastqc/{read_status}/{sample}_fastqc-data-{read_status}.txt"
    params:
        fname = lambda wc: wc.sample + "_" + wc.read_status,
        adapter = config["cutadapt"]["adapter"]
    wildcard_constraints:
        read_status="aligned_noPCRdup|unaligned"
    threads : config["threads"]
    log: "logs/fastqc/fastqc_{read_status}-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.read_status}) &> {log}
        (bedtools bamtofastq -fq qual_ctrl/fastqc/{wildcards.read_status}/{params.fname}.fastq -i {input}) &>> {log}
        (fastqc --adapters <(echo -e "adapter\t{params.adapter}") --nogroup --noextract -t {threads} -o qual_ctrl/fastqc/{wildcards.read_status} qual_ctrl/fastqc/{wildcards.read_status}/{params.fname}.fastq) &>> {log}
        (unzip -p qual_ctrl/fastqc/{wildcards.read_status}/{params.fname}_fastqc.zip {params.fname}_fastqc/fastqc_data.txt > {output}) &>> {log}
        (rm qual_ctrl/fastqc/{wildcards.read_status}/{params.fname}.fastq) &>> {log}
        """

fastqc_dict = {
        "per_base_qual":
        {   "title" : "Per base sequence quality",
            "fields": "base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus"
            } ,
        "per_tile_qual":
        {   "title" : "Per tile sequence quality",
            "fields": "tile\tbase\tmean\tsample\tstatus"
            },
        "per_seq_qual":
        {   "title" : "Per sequence quality scores",
            "fields": "quality\tcount\tsample\tstatus"
            },
        "per_base_seq_content":
        {   "title" : "Per base sequence content",
            "fields": "base\tg\ta\tt\tc\tsample\tstatus"
            },
        "per_seq_gc":
        {   "title" : "Per sequence GC content",
            "fields": "gc_content\tcount\tsample\tstatus"
            },
        "per_base_n":
        {   "title" : "Per base N content",
            "fields": "base\tn_count\tsample\tstatus"
            },
        "seq_length_dist":
        {   "title" : "Sequence Length Distribution",
            "fields": "length\tcount\tsample\tstatus"
            },
        "seq_duplication":
        {   "title" : "Total Deduplicated Percentage",
            "fields": "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus"
            },
        "adapter_content":
        {   "title" : "Adapter Content",
            "fields": "position\tpct\tsample\tstatus"
            }
        }

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}_fastqc-data-raw.txt", sample=SAMPLES),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}_fastqc-data-cleaned.txt", sample=SAMPLES),
        aligned_noPCRdup = expand("qual_ctrl/fastqc/aligned_noPCRdup/{sample}_fastqc-data-aligned_noPCRdup.txt", sample=SAMPLES),
        unaligned = expand("qual_ctrl/fastqc/unaligned/{sample}_fastqc-data-unaligned.txt", sample=SAMPLES),
    output:
        per_base_qual = 'qual_ctrl/fastqc/tss-seq-per_base_quality.tsv',
        per_tile_qual = 'qual_ctrl/fastqc/tss-seq-per_tile_quality.tsv',
        per_seq_qual =  'qual_ctrl/fastqc/tss-seq-per_sequence_quality.tsv',
        per_base_seq_content = 'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.tsv',
        per_seq_gc = 'qual_ctrl/fastqc/tss-seq-per_sequence_gc.tsv',
        per_base_n = 'qual_ctrl/fastqc/tss-seq-per_base_n.tsv',
        seq_length_dist = 'qual_ctrl/fastqc/tss-seq-sequence_length_distribution.tsv',
        seq_duplication = 'qual_ctrl/fastqc/tss-seq-sequence_duplication_levels.tsv',
        adapter_content = 'qual_ctrl/fastqc/tss-seq-adapter_content.tsv',
    run:
        shell("""rm -f {output}""")
        for fastqc_metric in output.keys():
            title = fastqc_dict[fastqc_metric]["title"]
            fields = fastqc_dict[fastqc_metric]["fields"]
            out_path = output[fastqc_metric]
            for read_status in input.keys():
                for sample_id, fastqc_data in zip(SAMPLES.keys(), input[read_status]):
                    shell("""awk 'BEGIN{{FS=OFS="\t"}} /{title}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{print $0, "{sample_id}", "{read_status}"}}' {fastqc_data} | tail -n +2 >> {out_path}""")
            shell("""sed -i "1i {fields}" {out_path}""")

rule plot_fastqc_summary:
    input:
        seq_len_dist = 'qual_ctrl/fastqc/tss-seq-sequence_length_distribution.tsv',
        per_tile = 'qual_ctrl/fastqc/tss-seq-per_tile_quality.tsv',
        per_base_qual = 'qual_ctrl/fastqc/tss-seq-per_base_quality.tsv',
        per_base_seq = 'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.tsv',
        per_base_n = 'qual_ctrl/fastqc/tss-seq-per_base_n.tsv',
        per_seq_gc = 'qual_ctrl/fastqc/tss-seq-per_sequence_gc.tsv',
        per_seq_qual = 'qual_ctrl/fastqc/tss-seq-per_sequence_quality.tsv',
        adapter_content = 'qual_ctrl/fastqc/tss-seq-adapter_content.tsv',
        seq_dup = 'qual_ctrl/fastqc/tss-seq-sequence_duplication_levels.tsv',
    output:
        seq_len_dist = 'qual_ctrl/fastqc/tss-seq-sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/tss-seq-per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/tss-seq-per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/tss-seq-per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/tss-seq-per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/tss-seq-per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/tss-seq-adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/tss-seq-sequence_duplication_levels.svg',
    script: "../scripts/fastqc_summary.R"

