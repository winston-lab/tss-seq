#!/usr/bin/env python

rule classify_peaks_genic:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-genic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.annotation} -wo -s | cut --complement -f11-13,15-17 > {output}
        """

rule classify_peaks_intragenic:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-intragenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.orfs} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}} $6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule classify_peaks_antisense:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        transcripts = config["genome"]["transcripts"]
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-antisense.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.transcripts} -wo -S | awk 'BEGIN{{FS=OFS="\t"}} $6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}}$6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule classify_peaks_convergent:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-convergent.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}} $6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}}$6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule classify_peaks_divergent:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-divergent.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -wo -s | awk 'BEGIN{{FS=OFS="\t"}}$6=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $12-$2+$10}}$6=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$13}}' > {output}
        """

rule classify_peaks_intergenic:
    input:
        peaks = "peakcalling/{group}/{group}-exp-idrpeaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "peakcalling/{group}/{group}-exp-idrpeaks-intergenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.transcripts} {input.orfs} -wa -v | bedtools intersect -a stdin -b {input.genic_anno} -wa -v -s | bedtools intersect -a stdin -b {input.annotation} -wa > {output}
        """

rule peakstats:
    input:
        expand("peakcalling/{group}/{group}-exp-idrpeaks-{category}.tsv", group=validgroups, category=CATEGORIES),
    output:
        table = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaknumbers.tsv",
        size = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaksizes-histogram.svg",
        violin_area = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaksizes-violin-equalarea.svg",
        violin_count = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peaksizes-violin-countscaled.svg",
        dist = "peakcalling/peakstats/{condition}-v-{control}/{condition}-v-{control}-peakdistances.svg"
    params:
        groups = lambda wc: [g for sublist in zip(controlgroups, conditiongroups) for g in sublist] if wc.condition=="all" else [wc.control, wc.condition]
    script:
        "scripts/peakstats.R"

rule classify_genic_diffexp_peaks:
    input:
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-genic-{direction}-summits.bed",
    log : "logs/classify_genic_diffexp_peaks/classify_genic_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b <(cut -f1-6 {input.annotation}) -wo -s | cut --complement -f21 | cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_score\ttranscript_name\ttranscript_score\ttranscript_strand")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.annotation} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_intragenic_diffexp_peaks:
    input:
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        orf_anno = config["genome"]["orf-annotation"],
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intragenic-{direction}-summits.bed",
    log : "logs/classify_intragenic_diffexp_peaks/classify_intragenic_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.orf_anno}) -wo -s | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=summit-$12}} $6=="-"{{$17=$13-(summit+1)}} {{print $0}}'  | cat <(paste <(head -n 1 {input.results}) <(echo -e "orf_chrom\torf_start\torf_end\torf_score\torf_name\torf_score\torf_strand\tatg_to_peak_dist")) - > {output.results}) &> {log}
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
    log : "logs/classify_antisense_diffexp_peaks/classify_antisense_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b <(cut -f1-6 {input.transcript_anno}) -wo -S | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=$13-(summit+1)}} $6=="-"{{$17=summit-$12}} {{print $0}}' | cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_score\ttranscript_name\ttranscript_strand\tsense_tss_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.transcript_anno} -u -S | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_convergent_diffexp_peaks:
    input:
        conv_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "convergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-convergent-{direction}-summits.bed",
    log : "logs/classify_convergent_diffexp_peaks/classify_convergent_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.conv_anno}) -wo -s | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=$13-(summit+1)}} $6=="-"{{$17=summit-$12}} {{print $0}}' | cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_score\ttranscript_name\ttranscript_strand\tsense_tss_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.conv_anno} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_divergent_diffexp_peaks:
    input:
        div_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "divergent-regions.bed",
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-divergent-{direction}-summits.bed",
    log : "logs/classify_divergent_diffexp_peaks/classify_divergent_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b <(cut -f1-6 {input.div_anno}) -wo -s | awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $6=="+"{{$17=summit-$12}} $6=="-"{{$17=$13-(summit+1)}} {{print $0}}' | cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_score\ttranscript_name\ttranscript_strand\tsense_tss_to_peak_dist")) - > {output.results}) &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.genic_anno} -v -s | bedtools intersect -a stdin -b {input.div_anno} -u -s | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

rule classify_intergenic_diffexp_peaks:
    input:
        intergenic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcript_anno = config["genome"]["transcripts"],
        orf_anno = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.narrowpeak",
        results = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{direction}.tsv",
    output:
        results = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intergenic-{direction}.tsv",
        narrowpeak = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intergenic-{direction}.narrowpeak",
        bed = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-intergenic-{direction}-summits.bed",
    log : "logs/classify_intergenic_diffexp_peaks/classify_intergenic_diffexp_peaks-{condition}-v-{control}_{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.results} -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | bedtools intersect -a stdin -b {input.intergenic_anno} -u | cat <(head -n 1 {input.results}) - > {output.results} &> {log}
        (bedtools intersect -a {input.narrowpeak} -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | bedtools intersect -a stdin -b {input.intergenic_anno} -u | tee {output.narrowpeak} | awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.bed}) &>> {log}
        """

#TODO: update these rules for all transcript classes
rule get_de_intragenic_frequency:
    input:
        orfs = config["genome"]["orf-annotation"],
        intrabed = "diff_exp/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}-results-{norm}-{direction}-intragenic.bed"
    output:
        "diff_exp/{condition}-v-{control}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
    shell: """
        bedtools intersect -a {input.orfs} -b {input.intrabed} -c -s > {output}
        """

rule plot_de_intragenic_frequency:
    input:
        "diff_exp/{condition}-v-{control}/{norm}/intragenic/intrafreq/{condition}-v-{control}-results-{norm}-{direction}-intrafreq.tsv"
    output:
        "diff_exp/{condition}-v-{control}/{norm}/intragenic/intrafreq/{condition}-v-{control}-intragenic-{norm}-{direction}-freqperORF.svg"
    log: "logs/plot_intragenic_frequency/plot_intragenic_frequency-{condition}-v-{control}-{norm}-{direction}.log"
    script: "scripts/intrafreq.R"

