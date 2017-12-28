library(tidyverse)
library(forcats)
library(gridExtra)
library(genefilter)
    
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                                    1))
        quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
        ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

main = function(groups, condition, control, table.out, size.out, violin.area.out, violin.count.out, dist.out){
    dflist = list()
    groups = unique(groups)
    for (i in 1:length(groups)){
        g = groups[i]
        dflist[[g]][['all']] = read_tsv(paste0('peakcalling/', groups[i], '-exp-idrpeaks.narrowPeak'), col_names=FALSE) %>%
            select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10) %>% distinct(.keep_all=TRUE)
        dflist[[g]][['genic']] = read_tsv(paste0('peakcalling/genic/', groups[i], '-exp-idrpeaks-genic.tsv'), col_names=FALSE) %>%
            select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11)
        dflist[[g]][['intragenic']] = read_tsv(paste0('peakcalling/intragenic/', groups[i], '-exp-idrpeaks-intragenic.tsv'), col_names=FALSE) %>%
            select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.ATG=X12) %>% distinct(chrom, start, end, strand, summit, .keep_all=TRUE) 
        dflist[[g]][['antisense']] = read_tsv(paste0('peakcalling/antisense/', groups[i], '-exp-idrpeaks-antisense.tsv'), col_names=FALSE) %>%
            select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.senseTSS=X12) %>% distinct(chrom, start, end, strand, summit, .keep_all=TRUE)
        dflist[[g]][['convergent']] = read_tsv(paste0('peakcalling/convergent/', groups[i], '-exp-idrpeaks-convergent.tsv'), col_names=FALSE) %>%
            select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.senseTSS=X12) %>% distinct(chrom, start, end, strand, summit, .keep_all=TRUE)
        dflist[[g]][['divergent']] = read_tsv(paste0('peakcalling/divergent/', groups[i], '-exp-idrpeaks-divergent.tsv'), col_names=FALSE) %>%
            select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.senseTSS=X12) %>% distinct(chrom, start, end, strand, summit, .keep_all=TRUE)
        dflist[[g]][['intergenic']] = read_tsv(paste0('peakcalling/intergenic/', groups[i], '-exp-idrpeaks-intergenic.tsv'), col_names=FALSE) %>%
            select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10) %>% distinct(.keep_all=TRUE)
    }
    
    ndf = tibble()
    for (i in 1:length(groups)){
        g = groups[i]
        dd = tibble(group=g,
                    all=nrow(dflist[[g]][['all']]),
                    genic=nrow(dflist[[g]][['genic']]),
                    intragenic=nrow(dflist[[g]][['intragenic']]),
                    antisense=nrow(dflist[[g]][['antisense']]),
                    convergent=nrow(dflist[[g]][['convergent']]),
                    divergent=nrow(dflist[[g]][['divergent']]),
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
        bind_rows(intergenicdf %>% transmute(size=end-start, group=group, type='intergenic')) %>% 
        mutate_at(vars(type,group), funs(fct_inorder(., ordered=TRUE)))
    
    sizedf %>% group_by(group, type) %>%
        summarise(n= n()) %>%
        spread(key=type, value=n) %>%
        ungroup() %>% 
        write_tsv(table.out)
    
    sizeannodf = sizedf %>% group_by(type, group, size) %>%
        mutate(freq = n()) %>%
        group_by(type) %>% mutate(y = max(freq)) %>%
        group_by(group, type) %>%
        summarise(size = .9*quantile(sizedf$size, .99), y = .9*unique(y), n = n()) 
    
    sizehist = ggplot() +
        geom_histogram(data = sizedf, aes(size), binwidth=1, center=0, fill="#08306b") +
        geom_text(data = sizeannodf, aes(x=size, y = y, label = n), hjust=1, vjust=1, size=4, fontface="bold") +
        facet_grid(type~group, scales="free_y", space="free_y", switch="y") +
        scale_x_continuous(limits = c(NA, quantile(sizedf$size, .99)), name="peak size (nt)") +
        scale_y_continuous(position="right", name=NULL, breaks=scales::pretty_breaks(n=4)) +
        ggtitle("TSS-seq peak sizes") +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=10, color="black"),
              strip.text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_text(angle=180, hjust=1),
              strip.background = element_blank(),
              panel.grid.major = element_line(color="grey80"),
              panel.grid.minor = element_line(color="grey80"))
    
    ggsave(size.out, plot = sizehist,
           width=length(groups)*8, height=16, units = "cm", limitsize=FALSE)
    
    # ggplot(data = sizedf, aes(x=size, y=fct_rev(type))) +
    #     geom_density_ridges(bandwidth=3) +
    #     scale_x_continuous(limits = c(NA, quantile(sizedf$size, .99)),
    #                        name="peak size (nt)",
    #                        expand=c(0.01, 0)) +
    #     scale_y_discrete(name=NULL, expand=c(0.02, 0)) +
    #     ggtitle("peak sizes") +
    #     facet_grid(.~group) +
    #     theme_ridges() +
    #     theme(axis.title.x = element_text(hjust=0.5, size=12, face="bold"),
    #           axis.text = element_text(size=12, face="bold"),
    #           panel.grid.major.x = element_line(color="grey40"),
    #           panel.grid.minor.x = element_line(color="grey40"),
    #           strip.background = element_blank(),
    #           strip.text = element_text(size=12, face="bold"))
        
    violin_base = ggplot() +
        scale_fill_brewer(palette="Set1", direction=-1) +
        scale_y_continuous(limits = c(0, quantile(sizedf$size, .99)), name="peak size (nt)") +
        theme_bw() +
        theme(legend.position="top",
              legend.title = element_blank(),
              legend.text = element_text(size=12, face="bold", color="black"),
              panel.grid.major.y = element_line(color="grey60"),
              panel.grid.minor.y = element_line(color="grey60"),
              panel.grid.major.x = element_blank(),
              axis.title.x = element_blank(),
              text = element_text(size=12, face="bold", color="black"),
              axis.text = element_text(size=12, face="bold", color="black"),
              axis.text.x = element_text(angle=30, hjust=.9))

    if(condition != "all"){
        violin_areascaled = violin_base +
            ggtitle(paste("peak sizes,", condition, "vs.", control)) +
            geom_split_violin(data = sizedf, aes(x=type, y=size, fill=group),
                              bw=3, scale="area")
        ggsave(violin.area.out, plot=violin_areascaled, width=14, height=10, units="cm")
        violin_countscaled = violin_base +
            ggtitle(paste("peak sizes,", condition, "vs.", control)) +
            geom_split_violin(data = sizedf %>% filter(type != "all"),
                              aes(x=type, y=size, fill=group), bw=3, scale="count")
        ggsave(violin.count.out, plot=violin_countscaled, width=14, height=10, units="cm")
    } else {
        violin_areascaled = violin_base +
            ggtitle("peak sizes") +
            geom_violin(data = sizedf, aes(x=type, y=size, fill=group),
                        bw=3, scale="area")
        ggsave(violin.area.out, plot=violin_areascaled, width=7*length(groups), height=10, units="cm")

        violin_countscaled = violin_base +
            ggtitle("peak sizes") +
            geom_violin(data = sizedf %>% filter(type != "all"),
                        aes(x=type, y=size, fill=group), bw=3, scale="count")
        ggsave(violin.count.out, plot=violin_countscaled, width=7*length(groups), height=10, units="cm")
    }
    
    make_annodf = function(df, colname, binsize){
        annodf = df %>% select(group, dist=colname) %>%
            mutate(x=quantile(dist, .995)) %>% group_by(group) %>%
            summarize(median=median(dist), max=max(dist), min=min(dist),
                      x=first(x), n=n(), mean_of_shorth=round(shorth(dist, tie.action="min"))) %>%
            mutate(y=.98*max(hist(df[[colname]], plot=FALSE,
                                  breaks=seq(min(min)-binsize, max(max)+binsize, binsize))$counts))
        return(annodf)
    }
    
    intra_annodf = make_annodf(intragenicdf, colname="peak.dist.to.ATG", binsize=25)
    as_annodf = make_annodf(antisensedf, colname="peak.dist.to.senseTSS", binsize=25) 
    conv_annodf = make_annodf(convergentdf, colname="peak.dist.to.senseTSS", binsize=10)
    div_annodf = make_annodf(divergentdf, colname="peak.dist.to.senseTSS", binsize=10)
    
    txn_dist_plot = function(df, annodf, colname, binsize){
        plot_base = ggplot() +
            geom_histogram(data = df, aes_string(colname), binwidth=binsize, boundary=0, fill="#3f007d") +
            geom_text(data = annodf, aes(x=x, y=y,
                                         label=paste0("n= ", n,
                                                      "\nmean of shorth= ", mean_of_shorth,
                                                      "\nmedian= ",median,
                                                      "\nmax= ", max)),
                      hjust=1, vjust=1, size=3, fontface="bold") +
            scale_y_continuous(name=NULL) +
            facet_grid(~group) +
            theme_light() +
            theme(text = element_text(size=12, color="black", face="bold"),
                  axis.text = element_text(size=10, color="black"),
                  strip.text = element_text(size=12, color="black", face="bold"),
                  strip.text.y = element_text(angle=180, hjust=1),
                  strip.background = element_blank(),
                  panel.grid.major = element_line(color="grey80"),
                  panel.grid.minor = element_line(color="grey80"),
                  plot.margin = unit(c(0.75,1.25,0.75,1.25), "cm"))
        return(plot_base)
    }
    
    intra_dist_plot = txn_dist_plot(intragenicdf, intra_annodf, "peak.dist.to.ATG", 25) +
        ggtitle("intragenic transcripts") +
        scale_x_continuous(limits = c(NA, quantile(intragenicdf$peak.dist.to.ATG, .995)),
                           name="distance: ATG to intragenic TSS (nt)",
                           minor_breaks = scales::pretty_breaks(n=10))
    
    as_dist_plot = txn_dist_plot(antisensedf, as_annodf, "peak.dist.to.senseTSS", 25) +
        ggtitle("all antisense transcripts") +
        scale_x_continuous(limits = c(NA, quantile(antisensedf$peak.dist.to.senseTSS, .995)),
                           name="distance: antisense TSS to sense TSS (bp)",
                           minor_breaks = scales::pretty_breaks(n=10))
    
    conv_dist_plot = txn_dist_plot(convergentdf, conv_annodf, "peak.dist.to.senseTSS", 10) +
        ggtitle("convergent antisense transcripts") +
        scale_x_continuous(limits = c(NA, quantile(convergentdf$peak.dist.to.senseTSS, .995)),
                           name="distance: convergent TSS to sense TSS (bp)",
                           minor_breaks = scales::pretty_breaks(n=10))
    
    div_dist_plot = txn_dist_plot(divergentdf, div_annodf, "peak.dist.to.senseTSS", 10) +
        ggtitle("divergent antisense transcripts") +
        scale_x_continuous(limits = c(NA, quantile(divergentdf$peak.dist.to.senseTSS, .995)),
                           name="distance: divergent TSS to sense TSS (bp)",
                           minor_breaks = scales::pretty_breaks(n=10))
    
    p = arrangeGrob(intra_dist_plot, as_dist_plot, conv_dist_plot, div_dist_plot, ncol=1)
    ggsave(dist.out, plot = p, width=length(groups)*10, height=26, units="cm")
}

main(groups = snakemake@params[["groups"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     table.out = snakemake@output[["table"]],
     size.out = snakemake@output[["size"]],
     violin.area.out = snakemake@output[["violin_area"]],
     violin.count.out = snakemake@output[["violin_count"]],
     dist.out = snakemake@output[["dist"]])
