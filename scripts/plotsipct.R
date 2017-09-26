library(tidyverse)

main = function(intable, barpath, boxpath, outstats){
    df = read_table2(intable, col_names=c('sample', 'group',
                                          'total','exp', 'si')) %>%
            mutate(sipct = si/total*100) %>%
            group_by(group) %>%
            mutate(outlier= ifelse(sipct>2.5*quantile(sipct,.75)-1.5*quantile(sipct,.25) |
                                    sipct< -2.5*quantile(sipct,.25)-1.5*quantile(sipct,.75),
                                    TRUE, FALSE)) %>%
            ungroup()
    
    nsamples = nrow(df)
    ngroups = df$group %>% unique %>% length()
    
    barplot = ggplot(data = df, aes(x=sample, fill=group, y = sipct)) +
                geom_col() +
                geom_text(aes(label=round(sipct, 1)),
                          size=5,
                          position=position_stack(vjust=0.9)) +
                scale_fill_brewer(palette='Set1') +
                ylab("% spike-in") +
                theme_minimal() +
                theme(axis.text = element_text(size=12, color="black"),
                      axis.text.x = element_text(angle=30, hjust=0.8),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=12, face="bold"),
                      legend.text = element_text(size=12),
                      legend.title = element_blank())
    
    boxplot = ggplot(data = df, aes(x=group, y=sipct, fill=group)) +
                geom_boxplot(outlier.size=0, outlier.stroke=0) +
                geom_point() +
                scale_fill_brewer(palette='Set1') +
                ylab("% spike-in") +
                theme_minimal() +
                theme(axis.text = element_text(size=12, color="black"),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=12),
                      legend.position = "none")
    
    ggsave(barpath, plot = barplot,
           width = 2+1.6*nsamples, height = 8, units = "cm")
    
    ggsave(boxpath, plot = boxplot,
           width = 2+2.5*ngroups, height = 8, units = "cm")
    
    outstats = df %>% add_count(group) %>% group_by(group) %>%
                mutate(median = median(sipct)) %>% ungroup() %>%
                filter(outlier==FALSE) %>% add_count(group) %>%
                group_by(group) %>%
                summarise(n = median(n), median = median(median),
                          n_no_outlier = median(nn),
                          mean_no_outlier = mean(sipct),
                          sd_no_outlier = sd(sipct)) %>%
                write.table(file = outstats, quote=FALSE, sep = "\t",
                            row.names=FALSE, col.names=TRUE)
}

df = main(intable = snakemake@input[[1]],
     barpath = snakemake@output[["barplot"]],
     boxpath = snakemake@output[["boxplot"]],
     outstats = snakemake@output[["stats"]])
