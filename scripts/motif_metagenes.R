library(tidyverse)
library(forcats)

main = function(in_table, condition, control, out_path){
    df = read_tsv(in_table,
                   col_names=c("direction", "factor", "txn_type",
                               "index", "position", "counts")) %>% 
        mutate_at(vars(txn_type), funs(fct_inorder(., ordered=TRUE))) %>%
        group_by(direction, factor, txn_type, position) %>%
        summarise(mean = mean(counts),
                  sem = sd(counts)/sqrt(n())) %>%
        ungroup() 
    
    nmotifs = df$factor %>% unique() %>% length()
    nclasses = df$txn_type %>% levels() %>% length()
    
    colors = c("up"="#b2182b", "unchanged"="grey30", "down"="#2166ac") 
    
    plot = ggplot(data = df, aes(x=position*1000, y=mean,
                                 ymin=pmax(0, mean-1.96*sem), ymax=mean+1.96*sem,
                                 color=fct_rev(direction),
                                 fill=fct_rev(direction))) +
        geom_vline(xintercept=0, color="grey50") +
        geom_ribbon(alpha=0.3, size=0) +
        geom_line(alpha=0.5) +
        facet_grid(factor~txn_type, scale="free_y", switch="y") +
        scale_color_manual(values=colors) +
        scale_fill_manual(values=colors) +
        scale_x_continuous(expand=c(0,0), name="position relative to peak start",
                           breaks=scales::pretty_breaks(n=3)) +
        ylab("mean motif frequency") +
        ggtitle("motif frequency at TSS-seq peaks",
                subtitle=paste(condition, "vs.", control)) +
        theme_light() +
        theme(strip.text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_text(angle=180, hjust=1),
              strip.background = element_blank(),
              strip.placement="outside",
              plot.title = element_text(size=12, color="black", face="bold"),
              plot.subtitle = element_text(size=12, color="black"),
              axis.title = element_text(size=12, color="black", face="bold"),
              legend.title = element_blank(),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size=12, color="black", face="bold"))
        
    ggsave(out_path, plot=plot, width=7+3*nclasses, height=4+1.5*nmotifs, units="cm",
           limitsize=FALSE)
}

main(in_table= snakemake@input[[1]],
     condition= snakemake@wildcards[["condition"]],
     control= snakemake@wildcards[["control"]],
     out_path= snakemake@output[[1]])
