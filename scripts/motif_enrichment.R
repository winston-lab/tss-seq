library(tidyverse)
library(broom)
library(ggrepel)
library(ggpmisc)

get_motif_counts = function(df){
    df %>% group_by(motif_id, motif_alt_id, tss_peak_id) %>% slice(1) %>% 
        group_by(motif_id, motif_alt_id) %>% count()
}

main = function(pval, alpha, condition, control, txn_type, direction,
                in_fimo_pos, in_fimo_neg, in_pos_total, in_neg_total,
                out_path, out_plot){
    fimo_positive = in_fimo_pos %>% read_tsv()
    pos_total = n_distinct(fimo_positive$tss_peak_id)
    fimo_positive = fimo_positive %>% filter(motif_logpadj>-log10(pval)) %>%
        get_motif_counts()
    
    fimo_negative = in_fimo_neg %>% read_tsv()
    neg_total = n_distinct(fimo_negative$tss_peak_id)
    fimo_negative = fimo_negative %>% filter(motif_logpadj>-log10(pval)) %>%
        get_motif_counts()
    
    df = full_join(fimo_positive, fimo_negative, by=c("motif_id", "motif_alt_id"),
                   suffix = c("_pos", "_neg")) %>%
        rename(pos_withmotif=n_pos, neg_withmotif=n_neg) %>% 
        mutate(pos_nomotif=pos_total-pos_withmotif,
               neg_nomotif=neg_total-neg_withmotif) %>%
        mutate_if(is.integer, funs(if_else(is.na(.), as.integer(0), .)))
    
    fisherdf = df %>% do(fisher.test(matrix(c(.$pos_withmotif, .$pos_nomotif,
                                              .$neg_withmotif, .$neg_nomotif),2,2),
                                     alternative="two.sided") %>% tidy()) %>% 
        rename(odds_ratio=estimate) %>%
        select(-c(method,alternative))
    
    df = df %>% left_join(fisherdf, by=c("motif_id", "motif_alt_id")) 
    df$fdr = p.adjust(df$p.value, method="BH")
    df = df %>% mutate_at(vars(odds_ratio), funs(log2(.))) %>% 
        arrange(fdr, p.value, desc(odds_ratio), desc(pos_withmotif)) %>% 
        select(motif_id, motif_alt_id, fdr, log2_odds_ratio=odds_ratio,
               conf_low=conf.low, conf_high=conf.high,
               pos_withmotif, pos_nomotif, neg_withmotif, neg_nomotif) %>% 
        mutate_if(is_double, funs(signif(., digits=3))) %>%
        write_tsv(out_path) %>%
        mutate(label=if_else(is.na(motif_alt_id), motif_id, motif_alt_id))
               
    plot = ggplot() +
        geom_hline(yintercept=-log10(alpha), linetype="dashed") +
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
    if(nrow(df %>% filter(fdr<alpha))>10){
        plot = plot +
            stat_dens2d_labels(data = df %>% filter(fdr<alpha),
                               aes(x=log2_odds_ratio, y=-log10(fdr), label=label),
                               geom="text_repel", keep.number=25, size=4)
    } else{
        plot = plot +
            geom_text_repel(data = df %>% filter(fdr<alpha),
                            aes(x=log2_odds_ratio, y=-log10(fdr), label=label),
                            size=4)
    }
    ggsave(out_plot, plot=plot, width=14, height=12, units="cm")
}
        

main(pval= snakemake@params[["pval_cutoff"]],
     alpha= snakemake@params[["alpha"]],
     condition= snakemake@wildcards[["condition"]],
     control= snakemake@wildcards[["control"]],
     txn_type = snakemake@wildcards[["category"]],
     direction= snakemake@params[["direction"]],
     in_fimo_pos = snakemake@input[["fimo_pos"]],
     in_fimo_neg = snakemake@input[["fimo_neg"]],
     out_path = snakemake@output[["tsv"]],
     out_plot = snakemake@output[["plot"]])
