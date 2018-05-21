library(tidyverse)
library(magrittr)
library(broom)
library(ggrepel)
library(ggpmisc)

#if motif is found more than once for a region,
#count only one
get_motif_counts = function(path, alpha, type="unchanged"){
    read_tsv(path, col_types = if(type=="unchanged"){'ciicdcdddiciicdccc'} else {'ciiiicciicdccc'}) %>%
        mutate(total_count = n_distinct(region_id)) %>%
        filter(motif_logpval > -log10(alpha)) %>%
        group_by(motif_id, motif_alt_id, region_id) %>%
        slice(1) %>%
        group_by(motif_id, motif_alt_id) %>%
        summarise(total_count=first(total_count), n=n()) %>%
        return()
}

main = function(fimo_pval, fdr_cutoff, condition, control, txn_type, direction, background_type,
                in_fimo_pos, in_fimo_neg, in_pos_total, in_neg_total,
                out_path, out_plot){
    df = get_motif_counts(in_fimo_pos, alpha=fimo_pval) %>%
        full_join(get_motif_counts(in_fimo_neg, alpha=fimo_pval, type=background_type),
                  by=c("motif_id", "motif_alt_id"),
                  suffix = c("_pos", "_neg")) %>%
        rename(pos_withmotif = n_pos, neg_withmotif = n_neg) %>%
        mutate(pos_nomotif = total_count_pos-pos_withmotif,
               neg_nomotif = total_count_neg-neg_withmotif) %>%
        mutate_if(is.integer, funs(if_else(is.na(.), as.integer(0), .))) %>%
        group_by(motif_id, motif_alt_id, pos_withmotif, neg_withmotif, pos_nomotif, neg_nomotif) %>%
        do(fisher.test(matrix(c(.$pos_withmotif, .$pos_nomotif,
                                .$neg_withmotif, .$neg_nomotif),2,2),
                       alternative="two.sided") %>% tidy()) %>%
        ungroup() %>%
        mutate(fdr = p.adjust(p.value, method="BH")) %>%
        mutate_at(vars(estimate, conf.low, conf.high), funs(log2(.))) %>%
        arrange(fdr, p.value, desc(estimate), desc(pos_withmotif)) %>%
        select(motif_id, motif_alt_id, fdr, log2_odds_ratio=estimate,
               conf_low=conf.low, conf_high=conf.high,
               pos_withmotif, pos_nomotif, neg_withmotif, neg_nomotif) %>%
        mutate_if(is_double, funs(signif(., digits=3))) %>%
        write_tsv(out_path) %>%
        mutate(label=if_else(is.na(motif_alt_id), motif_id, motif_alt_id))

    plot = ggplot() +
        geom_vline(xintercept=0, color="grey65") +
        geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed") +
        geom_point(data=df, aes(x=log2_odds_ratio, y=-log10(fdr)),
                   shape=16, size=1, alpha=0.8, stroke=0) +
        xlab(expression(bold(paste(log[2], " odds-ratio")))) +
        ylab(expression(bold(paste(-log[10], " FDR")))) +
        ggtitle(paste0("motif enrichment upstream of ", txn_type, " TSSs\n",
                       direction, " in ", condition, " vs. ", control),
                subtitle="Fisher's exact test (two-tailed)") +
        theme_light() +
        theme(text = element_text(size=12, face="bold", color="black"),
              axis.text = element_text(size=10, color="black"),
              plot.subtitle = element_text(face="plain"))
    if(nrow(df %>% filter(fdr < fdr_cutoff)) > 10){
        plot = plot +
            stat_dens2d_labels(data = df %>% filter(fdr < fdr_cutoff),
                               aes(x=log2_odds_ratio, y=-log10(fdr), label=label),
                               geom="text_repel", keep.number=25, size=4)
    } else{
        plot = plot +
            geom_text_repel(data = df %>% filter(fdr < fdr_cutoff),
                            aes(x=log2_odds_ratio, y=-log10(fdr), label=label),
                            size=4)
    }
    ggsave(out_plot, plot=plot, width=14, height=12, units="cm")
}

main(fimo_pval = snakemake@params[["fimo_pval"]],
     fdr_cutoff = snakemake@params[["fdr"]],
     condition= snakemake@wildcards[["condition"]],
     control= snakemake@wildcards[["control"]],
     txn_type = snakemake@wildcards[["category"]],
     direction= snakemake@params[["direction"]],
     background_type = snakemake@wildcards[["negative"]],
     in_fimo_pos = snakemake@input[["fimo_pos"]],
     in_fimo_neg = snakemake@input[["fimo_neg"]],
     out_path = snakemake@output[["tsv"]],
     out_plot = snakemake@output[["plot"]])

