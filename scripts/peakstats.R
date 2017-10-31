library(tidyverse)
library(forcats)
#library(modeest)
library(gridExtra)
    
main = function(groups, table.out, size.out, dist.out){
    dflist = list()
    for (i in 1:length(groups)){
        g = groups[i]
        dflist[[g]][['all']] = read_tsv(paste0('peakcalling/', groups[i], '-exp-idrpeaks.narrowPeak'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10)
        dflist[[g]][['genic']] = read_tsv(paste0('peakcalling/genic/', groups[i], '-exp-idrpeaks-genic.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11)
        dflist[[g]][['intragenic']] = read_tsv(paste0('peakcalling/intragenic/', groups[i], '-exp-idrpeaks-intragenic.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.ATG=X12)
        dflist[[g]][['antisense']] = read_tsv(paste0('peakcalling/antisense/', groups[i], '-exp-idrpeaks-antisense.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.senseTSS=X12)
        dflist[[g]][['convergent']] = read_tsv(paste0('peakcalling/convergent/', groups[i], '-exp-idrpeaks-convergent.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.senseTSS=X12)
        dflist[[g]][['divergent']] = read_tsv(paste0('peakcalling/divergent/', groups[i], '-exp-idrpeaks-divergent.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.senseTSS=X12)
        dflist[[g]][['intergenic']] = read_tsv(paste0('peakcalling/intergenic/', groups[i], '-exp-idrpeaks-intergenic.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10)
    }
    
    ndf = tibble()
    for (i in 1:length(groups)){
        g = groups[i]
        dd = tibble(group=g,
                    all=nrow(dflist[[g]][['all']]), genic=nrow(dflist[[g]][['genic']]),
                    intragenic=nrow(dflist[[g]][['intragenic']]), antisense=nrow(dflist[[g]][['antisense']]),
                    convergent=nrow(dflist[[g]][['convergent']]), divergent=nrow(dflist[[g]][['divergent']]),
                    intergenic=nrow(dflist[[g]][['intergenic']]))
        ndf = ndf %>% bind_rows(dd)
    }
    
    binddfs = function(type){
        df = tibble()
        for (i in 1:length(groups)){
            g = groups[i]
            dd = dflist[[g]][[type]] %>% mutate(group=g)
            df = df %>% bind_rows(dd)
        }
        df$group = fct_inorder(df$group, ordered=TRUE)
        return(df)
    }
    
    alldf = binddfs('all')
    genicdf = binddfs('genic')
    intragenicdf = binddfs('intragenic')
    antisensedf = binddfs('antisense')
    convergentdf = binddfs('convergent')
    divergentdf = binddfs('divergent')
    intergenicdf = binddfs('intergenic')
    
    sizedf = alldf %>% transmute(size=end-start, group=group, type='all') %>%
                bind_rows(genicdf %>% transmute(size=end-start, group=group, type='genic')) %>% 
                bind_rows(intragenicdf %>% transmute(size=end-start, group=group, type='intragenic')) %>% 
                bind_rows(antisensedf %>% transmute(size=end-start, group=group, type='antisense')) %>% 
                bind_rows(convergentdf %>% transmute(size=end-start, group=group, type='convergent')) %>% 
                bind_rows(divergentdf %>% transmute(size=end-start, group=group, type='divergent')) %>% 
                bind_rows(intergenicdf %>% transmute(size=end-start, group=group, type='intergenic'))
    sizedf$type = fct_inorder(sizedf$type, ordered=TRUE)
    sizedf$group = fct_inorder(sizedf$group, ordered=TRUE)
    
    sizedf %>% group_by(group, type) %>% summarise(n= n()) %>% spread(key=type, value=n) %>% ungroup() %>% 
        write_tsv(table.out)
    
    sizeannodf = sizedf %>% group_by(type, group, size) %>% mutate(freq = n()) %>%
                    group_by(type) %>% mutate(y = max(freq)) %>%
                    group_by(group, type) %>%
                    summarise(size = .7*quantile(sizedf$size, .9993), y = .5*unique(y), n = n()) 
    
    sizehist = ggplot() +
                geom_histogram(data = sizedf, aes(size), binwidth=1, fill="#08306b") +
                geom_text(data = sizeannodf, aes(x=size, y = y, label = n), hjust=1, size=4, fontface="bold") +
                facet_grid(type~group, scales="free_y", space="free_y", switch="y") +
                scale_x_continuous(limits = c(NA, quantile(sizedf$size, .9993)), name="peak size (nt)") +
                scale_y_continuous(position="right", name=NULL, breaks=scales::pretty_breaks(n=4)) +
                theme_light() +
                theme(text = element_text(size=12, color="black", face="bold"),
                      axis.text = element_text(size=10, color="black"),
                      strip.text = element_text(size=12, color="black", face="bold"),
                      strip.text.y = element_text(angle=180, hjust=1),
                      strip.background = element_blank(),
                      panel.grid.major = element_line(color="grey80"),
                      panel.grid.minor = element_line(color="grey80"))
    
    ggsave(size.out, plot = sizehist, width = length(groups)*7, height = 16, units = "cm", limitsize=FALSE)
    
    asannodf = antisensedf %>% select(group, dist = peak.dist.to.senseTSS) %>% mutate(x = .95*quantile(dist, .995)) %>% group_by(group) %>%
                #summarize(mode = mlv(dist, bw=100, method='parzen', kernel="gaussian")[['M']][[1]], median = median(dist), max=max(dist), x = first(x), n=n())
                summarize(median = median(dist), max=max(dist), x = first(x), n=n())
    
    asdistplot = ggplot() +
                    geom_histogram(data = antisensedf, aes(peak.dist.to.senseTSS), binwidth=25, fill="#08306b") +
                    geom_text(data = asannodf,
                              aes(x=x, label = paste0("n= ", n, "\nmedian= ",median, "\nmax= ", max)), y=60,
                              hjust=1, size=3, fontface="bold") +
                    scale_x_continuous(limits = c(NA, quantile(antisensedf$peak.dist.to.senseTSS, .995)),
                                       name="distance from antisense TSS to sense TSS (bp)",
                                       minor_breaks = scales::pretty_breaks(n=10)) +
                    scale_y_continuous(name=NULL) +
                    ggtitle("all antisense transcripts") +
                    facet_grid(~group) +
                    theme_light() +
                    theme(text = element_text(size=12, color="black", face="bold"),
                          axis.text = element_text(size=10, color="black"),
                          strip.text = element_text(size=12, color="black", face="bold"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.background = element_blank(),
                          panel.grid.major = element_line(color="grey80"),
                          panel.grid.minor = element_line(color="grey80"))
    
    convannodf = convergentdf %>% select(group, dist = peak.dist.to.senseTSS) %>% mutate(x = .95*quantile(dist, .995)) %>% group_by(group) %>%
                #summarize(mode = mlv(dist, bw=100, method='parzen', kernel="gaussian")[['M']][[1]], median = median(dist), max=max(dist), x = first(x), n=n())
                summarize(median = median(dist), max=max(dist), x = first(x), n=n())
    
    convdistplot = ggplot() +
                    geom_histogram(data = convergentdf, aes(peak.dist.to.senseTSS), binwidth=10, fill="#08306b") +
                    geom_text(data = convannodf,
                              aes(x=x, label = paste0("n= ", n, "\nmedian= ",median, "\nmax= ", max)), y=35,
                              hjust=1, size=3, fontface="bold") +
                    scale_x_continuous(limits = c(NA, quantile(convergentdf$peak.dist.to.senseTSS, .995)),
                                       name="distance from convergent TSS to sense TSS (bp)",
                                       minor_breaks = scales::pretty_breaks(n=10)) +
                    scale_y_continuous(name=NULL) +
                    ggtitle("convergent antisense transcripts") +
                    facet_grid(~group) +
                    theme_light() +
                    theme(text = element_text(size=12, color="black", face="bold"),
                          axis.text = element_text(size=10, color="black"),
                          strip.text = element_text(size=12, color="black", face="bold"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.background = element_blank(),
                          panel.grid.major = element_line(color="grey80"),
                          panel.grid.minor = element_line(color="grey80"))
    
    divannodf = divergentdf %>% select(group, dist = peak.dist.to.senseTSS) %>% mutate(x = .95*quantile(dist, .995)) %>% group_by(group) %>%
                #summarize(mode = mlv(dist, bw=100, method='parzen', kernel="gaussian")[['M']][[1]], median = median(dist), max=max(dist), x = first(x), n=n())
                summarize(median = median(dist), max=max(dist), x = first(x), n=n())
    
    divdistplot = ggplot() +
                    geom_histogram(data = divergentdf, aes(peak.dist.to.senseTSS), binwidth=10, fill="#08306b") +
                    geom_text(data = divannodf,
                              aes(x=x, label = paste0("n= ", n, "\nmedian= ",median, "\nmax= ", max)), y=40,
                              hjust=1, size=3, fontface="bold") +
                    scale_x_continuous(limits = c(NA, quantile(divergentdf$peak.dist.to.senseTSS, .995)),
                                       name="distance from divergent TSS to sense TSS (bp)",
                                       minor_breaks = scales::pretty_breaks(n=10)) +
                    scale_y_continuous(name=NULL) +
                    ggtitle("divergent antisense transcripts") +
                    facet_grid(~group) +
                    theme_light() +
                    theme(text = element_text(size=12, color="black", face="bold"),
                          axis.text = element_text(size=10, color="black"),
                          strip.text = element_text(size=12, color="black", face="bold"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.background = element_blank(),
                          panel.grid.major = element_line(color="grey80"),
                          panel.grid.minor = element_line(color="grey80"))
    
    p = arrangeGrob(asdistplot, convdistplot, divdistplot, ncol=1)
    ggsave(dist.out, plot = p, width = length(groups)*7, height=24, units="cm")
}

main(groups = snakemake@params[["groups"]],
     table.out = snakemake@output[["table"]],
     size.out = snakemake@output[["size"]],
     dist.out = snakemake@output[["dist"]])
