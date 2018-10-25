#!/usr/bin/env python

localrules: classify_genic_peaks, classify_intragenic_peaks, classify_antisense_peaks,
    classify_convergent_peaks, classify_divergent_peaks, classify_intergenic_peaks,
    classify_genic_diffexp_peaks, classify_intragenic_diffexp_peaks, classify_antisense_diffexp_peaks,
    classify_convergent_diffexp_peaks, classify_divergent_diffexp_peaks, classify_intergenic_diffexp_peaks,
    plot_peak_stats

peak_fields = "peak_chrom\tpeak_start\tpeak_end\tpeak_name\tpeak_score\tpeak_strand\tpeak_enrichment\tpeak_logpval\tpeak_logqval\tpeak_summit"

rule classify_genic_peaks:
    input:
        annotation = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        peaks = "peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak",
    output:
        table = "peakcalling/{group}/genic/{group}-experimental-idrpeaks-genic.tsv",
        narrowpeak = "peakcalling/{group}/genic/{group}-experimental-idrpeaks-genic.narrowpeak",
        bed = "peakcalling/{group}/genic/{group}-experimental-idrpeaks-genic-summits.bed",
    log : "logs/classify_peaks/classify_genic_peaks-{group}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | cut --complement -f17 | cat <(echo -e "{peak_fields}\tgenic_chrom\tgenic_start\tgenic_end\tgenic_name\tgenic_score\tgenic_strand") - > {output.table}) &> {log}
        (bedtools intersect -a {input.peaks} -b {input.annotation} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_intragenic_peaks:
    input:
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        orf_anno = config["genome"]["orf-annotation"],
        peaks = "peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak",
    output:
        table = "peakcalling/{group}/intragenic/{group}-experimental-idrpeaks-intragenic.tsv",
        narrowpeak = "peakcalling/{group}/intragenic/{group}-experimental-idrpeaks-intragenic.narrowpeak",
        bed = "peakcalling/{group}/intragenic/{group}-experimental-idrpeaks-intragenic-summits.bed",
    log : "logs/classify_peaks/classify_intragenic_peaks-{group}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.orf_anno}) -wo -s | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=summit-$12}} $6=="-"{{$17=$13-(summit+1)}} {{print $0}}' | cat <(echo -e "{peak_fields}\torf_chrom\torf_start\torf_end\torf_name\torf_score\torf_strand\tatg_to_peak_dist") - > {output.table}) &> {log}
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.orf_anno}) -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_antisense_peaks:
    input:
        transcript_anno = config["genome"]["transcripts"],
        peaks = "peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak",
    output:
        table = "peakcalling/{group}/antisense/{group}-experimental-idrpeaks-antisense.tsv",
        narrowpeak = "peakcalling/{group}/antisense/{group}-experimental-idrpeaks-antisense.narrowpeak",
        bed = "peakcalling/{group}/antisense/{group}-experimental-idrpeaks-antisense-summits.bed",
    log : "logs/classify_peaks/classify_antisense_peaks-{group}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b <(cut -f1-6 {input.transcript_anno}) -wo -S | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=$13-(summit+1)}} $6=="-"{{$17=summit-$12}} {{print $0}}' | cat <(echo -e "{peak_fields}\ttranscript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_peak_dist") - > {output.table}) &> {log}
        (bedtools intersect -a {input.peaks} -b <(cut -f1-6 {input.transcript_anno}) -u -S | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

# 0.) exclude genic peaks
# 1.) intersect with convergent region annotation
# 2.) join to original transcript annotation
# 3.) add distance information
rule classify_convergent_peaks:
    input:
        transcript_anno = config["genome"]["transcripts"],
        conv_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed"),
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        peaks = "peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak",
    output:
        table = "peakcalling/{group}/convergent/{group}-experimental-idrpeaks-convergent.tsv",
        narrowpeak = "peakcalling/{group}/convergent/{group}-experimental-idrpeaks-convergent.narrowpeak",
        bed = "peakcalling/{group}/convergent/{group}-experimental-idrpeaks-convergent-summits.bed",
    log : "logs/classify_peaks/classify_convergent_peaks-{group}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | sort -k14,14 | join -1 14 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4,4 {input.transcript_anno}) | sort -k1,1 -k2,2n | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=$13-(summit+1)}} $6=="-"{{$17=summit-$12}} {{print $0}}' | cat <(echo -e "{peak_fields}\ttranscript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_peak_dist") - > {output.table}) &> {log}
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

# 0.) exclude genic peaks
# 1.) intersect with divergent region annotation
# 2.) join to original transcript annotation
# 3.) add distance information
rule classify_divergent_peaks:
    input:
        transcript_anno = config["genome"]["transcripts"],
        div_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed"),
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        peaks = "peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak",
    output:
        table = "peakcalling/{group}/divergent/{group}-experimental-idrpeaks-divergent.tsv",
        narrowpeak = "peakcalling/{group}/divergent/{group}-experimental-idrpeaks-divergent.narrowpeak",
        bed = "peakcalling/{group}/divergent/{group}-experimental-idrpeaks-divergent-summits.bed",
    log : "logs/classify_peaks/classify_divergent_peaks-{group}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -wo -s | sort -k14,14 | join -1 14 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4,4 {input.transcript_anno}) | sort -k1,1 -k2,2n | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=(summit+1)-$13}} $6=="-"{{$17=$12-summit}} {{print $0}}' | cat <(echo -e "{peak_fields}\ttranscript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_peak_dist") - > {output.table}) &> {log}
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_intergenic_peaks:
    input:
        intergenic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"),
        transcript_anno = config["genome"]["transcripts"],
        orf_anno = config["genome"]["orf-annotation"],
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        peaks = "peakcalling/{group}/{group}_experimental-idrpeaks.narrowPeak",
    output:
        table = "peakcalling/{group}/intergenic/{group}-experimental-idrpeaks-intergenic.tsv",
        narrowpeak = "peakcalling/{group}/intergenic/{group}-experimental-idrpeaks-intergenic.narrowpeak",
        bed = "peakcalling/{group}/intergenic/{group}-experimental-idrpeaks-intergenic-summits.bed",
    log : "logs/classify_peaks/classify_intergenic_peaks-{group}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | bedtools intersect -a stdin -b {input.intergenic_anno} -u | cat <(echo -e "{peak_fields}") - > {output.table}) &> {log}
        (bedtools intersect -a {input.peaks} -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | bedtools intersect -a stdin -b {input.intergenic_anno} -u | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule plot_peak_stats:
    input:
        genic = "peakcalling/{group}/genic/{group}-experimental-idrpeaks-genic.tsv",
        intragenic = "peakcalling/{group}/intragenic/{group}-experimental-idrpeaks-intragenic.tsv",
        antisense = "peakcalling/{group}/antisense/{group}-experimental-idrpeaks-antisense.tsv",
        convergent = "peakcalling/{group}/convergent/{group}-experimental-idrpeaks-convergent.tsv",
        divergent = "peakcalling/{group}/divergent/{group}-experimental-idrpeaks-divergent.tsv",
        intergenic = "peakcalling/{group}/intergenic/{group}-experimental-idrpeaks-intergenic.tsv",
    output:
        table = "peakcalling/{group}/{group}-experimental-peak-stats.tsv",
        sizes = "peakcalling/{group}/{group}-experimental-peak-sizes.svg",
        distances = "peakcalling/{group}/{group}-experimental-peak-distances.svg",
    conda: "../envs/tidyverse.yaml"
    script:
        "../scripts/peak_stats.R"

rule classify_genic_diffexp_peaks:
    input:
        annotation = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}-summits.bed",
    log : "logs/classify_diffexp_peaks/classify_genic_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (tail -n +2 {input.results} | paste - <(cut -f10 {input.narrowpeak}) | bedtools intersect -a stdin -b {input.annotation} -wo -s | cut --complement -f22 | cat <(paste <(head -n 1 {input.results}) <(echo -e "peak_summit\tgenic_chrom\tgenic_start\tgenic_end\tgenic_name\tgenic_score\tgenic_strand")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.annotation} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_intragenic_diffexp_peaks:
    input:
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        orf_anno = config["genome"]["orf-annotation"],
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}-summits.bed",
    log : "logs/classify_diffexp_peaks/classify_intragenic_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (tail -n +2 {input.results} | paste - <(cut -f10 {input.narrowpeak}) | bedtools intersect -a stdin -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.orf_anno}) -wo -s | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$15}} $6=="+"{{$22=summit-$17}} $6=="-"{{$22=$18-(summit+1)}} {{print $0}}'  | cat <(paste <(head -n 1 {input.results}) <(echo -e "peak_summit\torf_chrom\torf_start\torf_end\torf_name\torf_score\torf_strand\tatg_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orf_anno} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_antisense_diffexp_peaks:
    input:
        transcript_anno = config["genome"]["transcripts"],
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-antisense-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-antisense-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-antisense-{direction}-summits.bed",
    log : "logs/classify_diffexp_peaks/classify_antisense_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (tail -n +2 {input.results} | paste - <(cut -f10 {input.narrowpeak}) | bedtools intersect -a stdin -b <(cut -f1-6 {input.transcript_anno}) -wo -S | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$15}} $6=="+"{{$22=$18-(summit+1)}} $6=="-"{{$22=summit-$17}} {{print $0}}' | cat <(paste <(head -n 1 {input.results}) <(echo -e "peak_summit\ttranscript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.transcript_anno} -u -S | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

# 0.) add summit information from narrowpeak to TSV
# 1.) exclude genic peaks
# 2.) intersect with convergent annotation
# 3.) add distance information
rule classify_convergent_diffexp_peaks:
    input:
        transcript_anno = config["genome"]["transcripts"],
        conv_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed"),
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-{direction}-summits.bed",
    log : "logs/classify_diffexp_peaks/classify_convergent_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (tail -n +2 {input.results} | paste - <(cut -f10 {input.narrowpeak}) | bedtools intersect -a stdin -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | sort -k19,19 | join -1 19 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4,4 {input.transcript_anno}) | sort -k11,11nr | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$15}} $6=="+"{{$22=$18-(summit+1)}} $6=="-"{{$22=summit-$17}} {{print $0}}' | cat <(paste <(head -n 1 {input.results}) <(echo -e "peak_summit\ttranscript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_divergent_diffexp_peaks:
    input:
        transcript_anno = config["genome"]["transcripts"],
        div_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed"),
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-{direction}-summits.bed",
    log : "logs/classify_diffexp_peaks/classify_divergent_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (tail -n +2 {input.results} | paste - <(cut -f10 {input.narrowpeak}) | bedtools intersect -a stdin -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.div_anno}) -wo -s | sort -k19,19 | join -1 19 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4,4 {input.transcript_anno}) | sort -k11,11nr | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$15}} $6=="+"{{$22=(summit+1)-$18}} $6=="-"{{$22=$17-summit}} {{print $0}}' | cat <(paste <(head -n 1 {input.results}) <(echo -e "peak_summit\ttranscript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_intergenic_diffexp_peaks:
    input:
        intergenic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"),
        transcript_anno = config["genome"]["transcripts"],
        orf_anno = config["genome"]["orf-annotation"],
        genic_anno = build_annotations(os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"),
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intergenic-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intergenic-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intergenic-{direction}-summits.bed",
    log : "logs/classify_diffexp_peaks/classify_intergenic_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (tail -n +2 {input.results} | paste - <(cut -f10 {input.narrowpeak}) | bedtools intersect -a stdin -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | bedtools intersect -a stdin -b {input.intergenic_anno} -u | cat <(paste <(head -n 1 {input.results}) <(echo "peak_summit")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | bedtools intersect -a stdin -b {input.intergenic_anno} -u | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

##TODO: update these rules for all transcript classes
#rule get_de_intragenic_frequency:
#    input:
#        orfs = config["genome"]["orf-annotation"],
#        intrabed = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic.bed"
#    output:
#        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
#    shell: """
#        bedtools intersect -a {input.orfs} -b {input.intrabed} -c -s > {output}
#        """

#rule plot_de_intragenic_frequency:
#    input:
#        "diff_exp/{condition}-v-{control}/{norm}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
#    output:
#        "diff_exp/{condition}-v-{control}/{norm}/intragenic/intrafreq/{condition}-v-{control}-intragenic-{norm}-{direction}-freqperORF.svg"
#    log: "logs/plot_intragenic_frequency/plot_intragenic_frequency-{condition}-v-{control}-{norm}-{direction}.log"
#    script: "scripts/intrafreq.R"

