library(tidyverse)
library(broom)
library(ggrepel)
library(ggpmisc)
library(viridis)
library(gridExtra)


plotclass = function(genicdf, classdf, condition, control, ttype){
    df = classdf %>% inner_join(genicdf, by="feat_name", suffix = c(".class", ".genic")) %>%
            group_by(feat_name) %>% arrange(desc(logpadj.class, logpadj.genic)) %>% slice(1) %>% ungroup()
    lfc.fit = lm(df$log2FoldChange.genic~df$log2FoldChange.class) %>% tidy()
    
    lfc.plot = ggplot(data = df, aes(x=log2FoldChange.class, y=log2FoldChange.genic)) +
            geom_vline(xintercept=0) +
            geom_hline(yintercept=0) +
            geom_smooth(method="lm", color="red") +
            geom_point(shape=16, alpha=0.6, size=0.5) +
            xlab(substitute(bold(type~log[2]~frac(num,den)), list(type=ttype, num=condition, den=control))) +
            ylab(substitute(bold(genic~log[2]~frac(num,den)), list(num=condition, den=control))) +
            ggtitle(paste("TSS-seq: genic vs.", ttype),
                    subtitle=paste("slope =", signif(lfc.fit$estimate[2],3),
                                   "(p =", signif(lfc.fit$p.value[2],3), ")")) +
            #stat_bin_hex(geom="point", aes(color=log10(..count..)), binwidth=c(0.1,0.1)) +
            #scale_color_viridis(option="inferno")
            stat_dens2d_filter(geom="text_repel", aes(label=feat_name), keep.fraction = 0.02, size=1) +
            theme_bw() +
            theme(axis.text = element_text(size=10),
                  axis.title = element_text(size=12, face="bold"),
                  axis.title.y = element_text(angle=0, vjust=0.5),
                  plot.title = element_text(size=12, face="bold"),
                  plot.subtitle = element_text(size=10))
            
    padj.plot = ggplot(data = df, aes(x=logpadj.class, y=logpadj.genic)) +
            geom_point(shape=16, alpha=0.6, size=0.5) +
            #stat_bin_hex(geom="point", aes(color=log10(..count..)), binwidth=c(0.8,0.8)) +
            #scale_color_viridis(option="inferno") +
            stat_dens2d_filter(geom="text_repel", aes(label=feat_name), keep.fraction = 0.01, size=1) + 
            xlab(substitute(bold(type:~-log[10]~p[adj]), list(type=ttype))) +
            ylab(expression(bold(genic:~-log[10]~p[adj]))) + 
            ggtitle("", subtitle="") + 
            theme_bw() +
            theme(axis.text = element_text(size=10),
                  axis.title = element_text(size=12, face="bold"),
                  plot.title = element_text(size=12, face="bold"),
                  plot.subtitle = element_text(size=10))
        
    return(list(lfc.plot, padj.plot, df))
}

main = function(condition, control, norm, genicdf, intradf, antidf, convdf, divdf, fig.path, df.path){
    genic = read_tsv(genicdf) %>% 
                select(meanExpr, log2FoldChange, logpadj, condition, control, feat_name=transcript_name)
    intra = read_tsv(intradf) %>% 
                select(peak_name, meanExpr, log2FoldChange, logpadj, condition, control, feat_name=ORF_name, dist=dist_ATG_to_peak)
    anti = read_tsv(antidf) %>% 
                select(peak_name, meanExpr, log2FoldChange, logpadj, condition, control, feat_name=transcript_name, dist=dist_peak_to_senseTSS)
    conv = read_tsv(convdf) %>% 
                select(peak_name, meanExpr, log2FoldChange, logpadj, condition, control, feat_name=transcript_name, dist=dist_peak_to_senseTSS)
    div = read_tsv(divdf) %>% 
                select(peak_name, meanExpr, log2FoldChange, logpadj, condition, control, feat_name=transcript_name, dist=dist_peak_to_senseTSS)
    
    ttypes = c("intragenic", "antisense", "convergent", "divergent")
    dfs = list(intra, anti, conv, div)
    
    plots = list()
    
    for (i in seq_along(ttypes)){
        plt = plotclass(genic, dfs[[i]], condition, control, ttypes[[i]])
        plots[[2*i-1]] = plt[[1]] 
        plots[[2*i]] = plt[[2]]
        write_tsv(plt[[3]],
                  path=paste0(df.path, condition, "-v-", control,
                              "-", norm, "-genic-v-", ttypes[[i]], ".tsv"))
    }
    
    out = arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                 plots[[5]], plots[[6]], plots[[7]], plots[[8]],
                 ncol=2, widths=unit(c(12,8), "cm"))
    
    ggsave(fig.path, out, height=30, width=20, units="cm")
}
main(condition=snakemake@wildcards[["condition"]],
     control=snakemake@wildcards[["control"]],
     norm=snakemake@wildcards[["norm"]],
     genicdf=snakemake@input[["genic"]],
     intradf=snakemake@input[["intragenic"]],
     antidf=snakemake@input[["antisense"]],
     convdf=snakemake@input[["convergent"]],
     divdf=snakemake@input[["divergent"]],
     fig.path=snakemake@output[["figure"]],
     df.path=snakemake@params[["path"]])
