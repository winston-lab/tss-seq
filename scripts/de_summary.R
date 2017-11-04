library(tidyverse)
library(forcats)
library(ggrepel)
library(ggpmisc)

alpha = snakemake@params[["alpha"]]

import = function(path){
    read_tsv(path) %>%
    mutate(sig = if_else(logpadj > -log10(alpha), TRUE, FALSE)) %>% 
    return()
}

clean = function(df, labelcol, type){
    df %>% select(label = labelcol, meanExpr, log2FoldChange, logpadj, sig) %>% mutate(type=type) %>% return()
}

bin = function(df, type){
    df %>%
    transmute(change = if_else(sig, if_else(log2FoldChange > 0, "up", "down"), "unchanged")) %>%
    group_by(change) %>% count() %>% filter(!is.na(change)) %>% mutate(class=type) %>% return()
}

main = function(in.all, in.genic, in.intra, in.as, in.conv, in.div, in.inter, cond, ctrl, lfc, alpha, out.ma, out.volcano, out.summary){
    all = import(in.all) 
    genic = import(in.genic)
    intra = import(in.intra)
    as = import(in.as)
    conv = import(in.conv)
    div = import(in.div)
    inter = import(in.inter)
    
    cleandf = clean(all, 'peak_name', 'all') %>% bind_rows(clean(genic, 'transcript_name', 'genic')) %>% 
                bind_rows(clean(intra, 'ORF_name', 'intragenic')) %>% bind_rows(clean(as, 'transcript_name', 'antisense')) %>% 
                bind_rows(clean(conv, 'transcript_name', 'convergent')) %>% bind_rows(clean(div, 'transcript_name', 'divergent')) %>% 
                bind_rows(clean(inter, 'peak_name', 'intergenic'))
    cleandf$type = fct_inorder(cleandf$type, ordered=TRUE)
    minx = quantile(cleandf$meanExpr, .2)
    
    maplot = ggplot(data = cleandf, aes(x=meanExpr, y = log2FoldChange)) +
                geom_hline(yintercept = 0, linetype="dashed") +
                geom_point(aes(color=sig), alpha=0.4, stroke=0, size=1) +
                scale_color_manual(values = c("grey40", "red"), guide=FALSE) +
                stat_dens2d_filter(data = cleandf %>% filter(sig & meanExpr>minx & log2FoldChange > 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_y = .3,
                                   size=1.5) +
                stat_dens2d_filter(data = cleandf %>% filter(sig & meanExpr>minx & log2FoldChange < 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_y = -.3,
                                   size=1.5) + 
                scale_x_log10(name="mean expression level") +
                scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
                ylab(substitute(bold(log[bold(2)]~frac(cond,cont)), list(cond=cond, cont=ctrl))) +
                facet_wrap(~type) +
                ggtitle(paste("TSS-seq MA plots:", cond, "vs.", ctrl),
                        subtitle = bquote("DESeq2: |" ~ log[2] ~ "fold-change" ~ "| >" ~ .(lfc) ~ "@ FDR" ~ .(alpha)))  +
                theme_bw() +
                theme(text = element_text(size=12, face="bold"),
                      axis.text = element_text(size=10, color="black"),
                      axis.title.y = element_text(angle=0, vjust=0.5),
                      strip.background = element_blank(),
                      strip.text = element_text(size=12, face="bold", color="black"))
    
    ggsave(out.ma, maplot, height=16, width=22, units="cm")
    
    volcano = ggplot(data = cleandf, aes(x = log2FoldChange, y = logpadj))+
                geom_point(aes(color=sig), alpha=0.4, stroke=0, size=1) +
                scale_color_manual(values = c("grey40", "red"), guide=FALSE) +
                stat_dens2d_filter(data = cleandf %>% filter(sig & log2FoldChange < 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_x = -0.3,
                                   size=1.5) +
                stat_dens2d_filter(data = cleandf %>% filter(sig & log2FoldChange > 0),
                                   geom = "text_repel",
                                   aes(label = label),
                                   keep.number = 5,
                                   point.padding = unit(0.1, "lines"),
                                   box.padding = unit(0.1, "lines"),
                                   nudge_x = 0.3,
                                   size=1.5) +
                ylab(expression(bold(log[10] ~ p[adj]))) +
                xlab(substitute(bold(log[bold(2)]~frac(cond,cont)), list(cond=cond, cont=ctrl))) +
                facet_wrap(~type) +
                ggtitle(paste("TSS-seq volcano plots:", cond, "vs.", ctrl),
                        subtitle = bquote("DESeq2: |" ~ log[2] ~ "fold-change" ~ "| >" ~ .(lfc) ~ "@ FDR" ~ .(alpha)))  +
                theme_bw() +
                theme(text = element_text(size=12, face="bold"),
                      axis.text = element_text(size=10, color="black"),
                      axis.title.y = element_text(angle=0, vjust=0.5),
                      strip.background = element_blank(),
                      strip.text = element_text(size=12, face="bold", color="black"))
    
    ggsave(out.volcano, volcano, height=16, width=20, units="cm")
    
    countdf = bind_rows(bin(genic, 'genic')) %>% bind_rows(bin(intra, 'intragenic')) %>%
                bind_rows(bin(as, 'antisense')) %>% bind_rows(bin(conv, 'convergent')) %>%
                bind_rows(bin(div, 'divergent')) %>% bind_rows(bin(inter, 'intergenic')) %>%
                group_by(class) %>% mutate(ymax = cumsum(n), ymin = (cumsum(n)-n)) %>%
                mutate_at(vars(ymin, ymax), funs(./max(ymax)))
    
    countdf$class = fct_inorder(countdf$class, ordered=TRUE)
    countdf$change = fct_inorder(countdf$change, ordered=TRUE)
    
    csx = countdf %>% group_by(class) %>% summarise(classn = sum(n)) %>%
            mutate(xmax = cumsum(classn), xmin = cumsum(classn)-classn) %>% select(-classn)
    countdf = countdf %>% left_join(csx, by='class')
    
    gsummary = ggplot() +
                geom_rect(data=countdf, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=change),
                          color="white", alpha=0.9, size=1.5) +
                geom_text(data = countdf, aes(x=(xmin+xmax)/2, y=(ymin+ymax)/2, label=n),
                          size=4, color="black", fontface="bold") + 
                geom_text(data = (countdf %>% summarise(x = (max(xmax)+min(xmin))/2)),
                            aes(x=x, label=class), y=1.04, angle=15, hjust=.2, size=4, fontface="bold") +
                scale_fill_brewer(palette = "Set1") +
                scale_x_continuous(limits = c(NA, max(countdf$xmax)*1.04)) +
                scale_y_continuous(limits = c(NA, 1.1)) +
                guides(fill=guide_legend(reverse=TRUE)) +
                theme_void() +
                ggtitle(paste("TSS-seq differential expression:", cond, "vs.", ctrl),
                        subtitle = bquote("DESeq2: |" ~ log[2] ~ "fold-change" ~ "| >" ~ .(lfc) ~ "@ FDR" ~ .(alpha)))  +
                theme(legend.text = element_text(size=12, face="bold", color="black"),
                      legend.title = element_blank(),
                      legend.position = c(.94, .5),
                      legend.justification = "left",
                      legend.key.size = unit(1, "cm"),
                      plot.margin = unit(c(.5, 3, 0, 0), "cm"),
                      plot.title = element_text(size=12, face="bold", color="black", hjust=0.16),
                      plot.subtitle = element_text(size=10, hjust=0.08))
    
    ggsave(out.summary, gsummary, height=10, width=20, units="cm")
}

main(in.all = snakemake@input[["total"]],
     in.genic = snakemake@input[["genic"]],
     in.intra = snakemake@input[["intragenic"]],
     in.as = snakemake@input[["antisense"]],
     in.conv = snakemake@input[["convergent"]],
     in.div = snakemake@input[["divergent"]],
     in.inter = snakemake@input[["intergenic"]],
     cond = snakemake@wildcards[["condition"]],
     ctrl = snakemake@wildcards[["control"]],
     lfc = snakemake@params[["lfc"]],
     alpha = snakemake@params[["alpha"]],
     out.ma = snakemake@output[["maplot"]],
     out.volcano = snakemake@output[["volcano"]],
     out.summary = snakemake@output[["summary"]])
