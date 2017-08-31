library(tidyverse)

raw = read_table2(snakemake@input[[1]],
                  col_names = c('chrom', 'strand', 'peakstart', 'peakend', 'peakname',
                                'orfstart','orfend','orfname','lfc','sig','disttogenic'))

df = raw %>% group_by(orfname, orfstart, orfend) %>% summarise(nintra = n())

plot = ggplot(data = df, aes(nintra)) +
        geom_histogram(binwidth=1, color="black", fill="#377eb8") +
        xlab("number of intragenic TSS' per\nORF with >0 intragenic TSS") +
        ylab(label=NULL) +
        ggtitle(paste0(snakemake@wildcards[["condition"]], "-v-", snakemake@wildcards[["control"]])) +
        theme_minimal() +
        theme(axis.title = element_text(size=12, face="bold"),
              axis.text = element_text(size=12),
              plot.title = element_text(size=12, face="bold"))

ggsave(snakemake@output[[1]], plot = plot, width = 9, height=7, units="cm")
