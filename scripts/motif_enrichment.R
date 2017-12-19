library(tidyverse)
library(broom)

import_fimo = function(path, alpha){
    read_tsv(path, skip=1,
             col_names=c('motif_id','motif_alt_id','seq_name',
                         'start','end','strand','score', 'pval','qval')) %>% 
        filter(qval<alpha)
}

get_motif_counts = function(df){
    df %>% group_by(motif_id, motif_alt_id, seq_name) %>% slice(1) %>% 
        group_by(motif_id, motif_alt_id) %>% count()
}

main = function(alpha, in_fimo_pos, in_fimo_neg, in_pos_total, in_neg_total, out_path){
    fimo_positive = in_fimo_pos %>% import_fimo(alpha=alpha) %>% get_motif_counts()
    fimo_negative = in_fimo_neg %>% import_fimo(alpha=alpha) %>% get_motif_counts()
    pos_total = read_tsv(in_pos_total, col_names=FALSE) %>% nrow()
    neg_total = read_tsv(in_neg_total, col_names=FALSE) %>% nrow()
    
    df = full_join(fimo_positive, fimo_negative, by=c("motif_id", "motif_alt_id"),
                   suffix = c("_pos", "_neg")) %>%
        rename(pos_withmotif=n_pos, neg_withmotif=n_neg) %>% 
        mutate(pos_nomotif=pos_total-pos_withmotif,
               neg_nomotif=neg_total-neg_withmotif)
    df[is.na(df)]=0
    
    fisherdf = df %>% do(fisher.test(matrix(c(.$pos_withmotif, .$pos_nomotif,
                                              .$neg_withmotif, .$neg_nomotif),2,2),
                                     alternative="g") %>% tidy()) %>% 
        rename(odds_ratio=estimate) %>%
        select(-c(method,alternative))
    
    df = df %>% left_join(fisherdf, by=c("motif_id", "motif_alt_id"))
    df$fdr = p.adjust(df$p.value, method="BH")
    df = df %>% arrange(fdr, p.value, desc(odds_ratio), desc(pos_withmotif)) %>% 
        write_tsv(out_path)
}

main(alpha= snakemake@params[["alpha"]],
     in_fimo_pos = snakemake@input[["fimo_pos"]],
     in_fimo_neg = snakemake@input[["fimo_neg"]],
     in_pos_total = snakemake@input[["pos_total"]],
     in_neg_total = snakemake@input[["neg_total"]],
     out_path = snakemake@output[[1]])
