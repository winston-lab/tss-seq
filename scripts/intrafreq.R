library(tidyverse)
library(gridExtra)

raw = read_tsv(snakemake@input[[1]],
                  col_names = c('chrom', 'orfstart','orfend','orfname','score','strand','nintra'))

p1= ggplot(data = raw, aes(nintra)) +
        geom_histogram(binwidth=1, color="black", fill="#377eb8") +
        xlab("intragenic TSSs per ORF (all ORFs)") +
        ylab(label=NULL) +
        ggtitle(paste0(snakemake@wildcards[["condition"]], "-v-", snakemake@wildcards[["control"]])) +
        theme_minimal() +
        theme(axis.title = element_text(size=12, face="bold"),
              axis.text = element_text(size=12),
              plot.title = element_text(size=12, face="bold"))

p2= ggplot(data = raw %>% filter(nintra >0), aes(nintra)) +
        geom_histogram(binwidth=1, color="black", fill="#377eb8") +
        xlab("intragenic TSSs per ORF\n(ORFs with intragenic)") +
        ylab(label=NULL) +
        #ggtitle(paste0(snakemake@wildcards[["condition"]], "-v-", snakemake@wildcards[["control"]])) +
        theme_minimal() +
        theme(axis.title = element_text(size=12, face="bold"),
              axis.text = element_text(size=12),
              plot.title = element_text(size=12, face="bold"))

plots = arrangeGrob(p1, p2)

ggsave(snakemake@output[[1]], plot = plots, width = 9, height=12, units="cm")
