library(tidyverse)
library(forcats)
library(gridExtra)

main = function(intable, samplelist, controls, conditions, plotpath, statspath){
    df = read_table2(intable, col_names=c('sample', 'group',
                                          'total','exp', 'si')) %>%
            filter(sample %in% samplelist) %>%
            mutate(sipct = si/total*100) %>%
            group_by(group) %>%
            mutate(outlier= ifelse(sipct>2.5*quantile(sipct,.75)-1.5*quantile(sipct,.25) |
                                    sipct< -2.5*quantile(sipct,.25)-1.5*quantile(sipct,.75),
                                    TRUE, FALSE)) %>%
            ungroup()
    df$sample = fct_inorder(df$sample, ordered=TRUE)
    df$group = fct_inorder(df$group, ordered=TRUE)
    
    nsamples = nrow(df)
    ngroups = df$group %>% unique %>% length()
    
    barplot = ggplot(data = df, aes(x=sample, fill=group, y = sipct)) +
                geom_col() +
                geom_text(aes(label=round(sipct, 1)),
                          size=4,
                          position=position_stack(vjust=0.9)) +
                scale_fill_brewer(palette='Set1', direction=-1) +
                ylab("% spike-in") +
                theme_minimal() +
                theme(axis.text = element_text(size=8, color="black"),
                      axis.text.x = element_text(angle=30, hjust=0.8),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=8, face="bold"),
                      legend.position = "none")
    
    boxplot = ggplot(data = df, aes(x=group, y=sipct, fill=group)) +
                geom_boxplot(outlier.size=1.5, outlier.color="red") +
                geom_point(size=1) +
                scale_fill_brewer(palette='Set1', direction=-1) +
                ylab("% spike-in") +
                theme_minimal() +
                theme(axis.text = element_text(size=8, color="black"),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size=8, face="bold"),
                      legend.position = "none")
    
    outstats = df %>% add_count(group) %>% group_by(group) %>%
                mutate(median = median(sipct)) %>% ungroup() %>%
                filter(outlier==FALSE) %>% add_count(group) %>%
                group_by(group) %>%
                summarise(n = median(n), median = median(median),
                          n_no_outlier = median(nn),
                          mean_no_outlier = mean(sipct),
                          sd_no_outlier = sd(sipct)) 
    write_tsv(outstats, path = statspath, col_names=TRUE)
    
    ctrl = outstats %>% inner_join(data_frame(c=controls), by=c("group"="c"))
    cond = outstats %>% inner_join(data_frame(c=conditions), by=c("group"="c"))
    pct = bind_cols(ctrl, cond) %>%
            select(group, group1, mean_no_outlier, mean_no_outlier1) %>%
            rename(control=group, condition=group1,
                   pct.ctrl = mean_no_outlier, pct.cond = mean_no_outlier1) %>%
            mutate_at(c("pct.ctrl", "pct.cond"), funs(./100)) %>%
            mutate(relative.levels = pct.ctrl*(1-pct.cond)/(pct.cond*(1-pct.ctrl)))
    
    pctout = pct %>% select(control, condition, relative.levels) %>%
                mutate_at("relative.levels", funs(round(., digits=3)))
    pctdraw = tableGrob(pctout, rows=NULL, cols=c("control","condition","relative levels"),
                        ttheme_default(base_size=8))
    #set width
    wl = 1+nsamples
    wr = 1+2*ngroups
    th = 1+ngroups/4
    page = arrangeGrob(barplot, boxplot, pctdraw, ncol=2,
                        widths=unit(c(wl, wr), c("cm","cm")),
                        heights=unit(c(7,th,0,0), c("cm","cm","cm","cm")))
    
    ggsave(plotpath, page, width = wl+wr, height=7+th+.5, units = "cm")
}

df = main(intable = snakemake@input[[1]],
          samplelist = snakemake@params[["samplelist"]],
          controls = snakemake@params[["controls"]],
          conditions = snakemake@params[["conditions"]],
          plotpath = snakemake@output[["plot"]],
          statspath = snakemake@output[["stats"]])
