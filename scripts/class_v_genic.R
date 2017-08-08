library(tidyverse)
library(ggrepel)

main = function(in.class, in.pclass.up, in.pclass.dn, in.genic, out.scatter.text, out.scatter.nolabel, out.table){
    pclass.colnames = c('chrom','strand','peak.start','peak.end','peak.name','feat.start','feat.end','feat.name','peak.lfc','peak.sig')
    pclass.up = read_table2(in.pclass.up, col_names=pclass.colnames)
    pclass.dn = read_table2(in.pclass.dn, col_names=pclass.colnames) %>%
                    mutate_at("peak.lfc", funs(-.))
    pclass = pclass.up %>% bind_rows(pclass.dn) %>% select(feat.name, peak.lfc, peak.sig)
    
    genic = read_table2(in.genic, col_names=TRUE) %>% select(base, log2FoldChange, padj)
    
    df = left_join(pclass, genic, by=c("feat.name" = "base")) %>% mutate_at("padj", funs(-log10(.))) %>%
            rename(name=feat.name, genic.lfc=log2FoldChange, genic.sig = padj) %>% drop_na()
    
    write.table(df, file=out.table, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    
    scatter.text = ggplot(data = df, aes(x=peak.lfc, y=genic.lfc)) +
                        geom_abline(slope=1, intercept=0, linetype="dashed", alpha=0.5) +
                        geom_point(size=0.5, alpha=0.8) +
                        geom_text_repel(aes(size=genic.sig, label=name), segment.alpha=0.7, box.padding=unit(.1, "lines")) +
                        xlab(expression(paste(in.class, " ", log[2], "(fold-change)"))) +
                        ylab(expression(paste("genic ", log[2], "(fold-change)"))) +
                        labs(size=expression(paste(-log[10], "(genic pval)"))) +
                        theme_bw() +
                        theme(axis.title=element_text(size=12, face="bold"),
                              legend.title = element_text(size=12, face="bold"),
                              axis.text=element_text(size=12))
    
    ggsave(out.scatter.text, plot=scatter.text, height=12, width=24, units="cm")
    
    scatter.nolabel = ggplot(data = df, aes(x=peak.lfc, y=genic.lfc)) +
                        geom_abline(slope=1, intercept=0, linetype="dashed", alpha=0.5) +
                        geom_point(aes(size=genic.sig), alpha=0.8) +
                        xlab(expression(paste(in.class, " ", log[2], "(fold-change)"))) +
                        ylab(expression(paste("genic ", log[2], "(fold-change)"))) +
                        labs(size=expression(paste(-log[10], "(genic pval)"))) +
                        theme_bw() +
                        theme(axis.title=element_text(size=12, face="bold"),
                              legend.title = element_text(size=12, face="bold"),
                              axis.text=element_text(size=12))
    
    ggsave(out.scatter.nolabel, plot=scatter.nolabel, height=12, width=24, units="cm")
}

main(in.class = snakemake@wildcards[["type"]],
     in.pclass.up = snakemake@input[["pclass_up"]],
     in.pclass.dn = snakemake@input[["pclass_dn"]],
     in.genic = snakemake@input[["genic"]],
     out.scatter.text = snakemake@output[["scatter_text"]],
     out.scatter.nolabel = snakemake@output[["scatter_nolabel"]],
     out.table = snakemake@output[["table"]])
