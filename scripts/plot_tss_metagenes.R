#!/usr/bin/env Rscript
library(psych)
library(tidyverse)
library(forcats)
library(viridis)
library(ggthemes)

#stat_stepribbon by user 'felasa':
#https://groups.google.com/forum/?fromgroups=#!topic/ggplot2/9cFWHaH1CPs
stairstepn <- function( data, direction="hv", yvars="y" ) {
    direction <- match.arg( direction, c( "hv", "vh" ) )
    data <- as.data.frame( data )[ order( data$x ), ]
    n <- nrow( data )
    
    if ( direction == "vh" ) {
        xs <- rep( 1:n, each = 2 )[ -2 * n ]
        ys <- c( 1, rep( 2:n, each = 2 ) )
    } else {
        ys <- rep( 1:n, each = 2 )[ -2 * n ]
        xs <- c( 1, rep( 2:n, each = 2))
    }
    
    data.frame(
        x = data$x[ xs ]
        , data[ ys, yvars, drop=FALSE ]
        , data[ xs, setdiff( names( data ), c( "x", yvars ) ), drop=FALSE ]
    ) 
}

stat_stepribbon <- 
    function(mapping = NULL, data = NULL, geom = "ribbon", position = "identity", inherit.aes = TRUE) {
        ggplot2::layer(
            stat = Stepribbon, mapping = mapping, data = data, geom = geom, 
            position = position, inherit.aes = inherit.aes
        )
    }

Stepribbon <- 
    ggproto("stepribbon", Stat,
            compute_group = function(., data, scales, direction = "hv", yvars = c( "ymin", "ymax" ), ...) {
                stairstepn( data = data, direction = direction, yvars = yvars )
            },                        
            required_aes = c( "x", "ymin", "ymax" )
    )

theme_default = theme_light() +
    theme(text = element_text(size=12, color="black", face="bold"),
          strip.background = element_blank(),
          strip.text = element_text(size=12, color="black", face="bold"),
          strip.text.y = element_text(angle=-180),
          axis.text.x = element_text(size=12, color="black", face="bold"),
          axis.text.y = element_text(size=10, color="black", face="plain"),
          axis.title = element_text(size=10, face="plain"),
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=12, face="plain"))
    
main = function(intable, samplelist, strand, upstream, dnstream, trim_pct,
                refptlabel, scaled_length, endlabel, ylabel, meta_sample_out, 
                meta_sample_overlay_out, meta_heatmap_sample_out, meta_group_out,
                meta_group_overlay_out, meta_heatmap_group_out){
    meta = function(df){
        ggplot(data = df, aes(x=position, y=mean, ymin = mean-1.96*sem,
                                    ymax=mean+1.96*sem, alpha=0.8)) +
        geom_vline(xintercept = c(0, scaled_length/1000), color="grey70") +
        stat_stepribbon() +
        geom_step(alpha=1, color="#114477") +
        scale_alpha(guide=FALSE) +
        scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                           labels=c(refptlabel, "", endlabel),
                           name="scaled distance",
                           limits = c(-upstream/1000, (dnstream+scaled_length)/1000),
                           expand=c(0,0)) +
        ylab("normalized counts") +
        ggtitle(paste("mean", strand, "TSS-seq signal"),
                subtitle = paste(nindices, ylabel)) +
        theme_default
    }
    meta_overlay = function(df){
        plot = ggplot(data = df,
                      aes(x=position, y=mean, ymin = mean-1.96*sem,
                          ymax=mean+1.96*sem, alpha=0.8, fill=group)) +
        geom_vline(xintercept = c(0, scaled_length/1000), color="grey70") +
        scale_alpha(guide=FALSE) +
        scale_fill_ptol(name=NULL) +
        scale_color_ptol(name=NULL) +
        scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                           labels=c(refptlabel, "", endlabel),
                           name="scaled distance",
                           limits = c(-upstream/1000, (dnstream+scaled_length)/1000),
                           expand=c(0,0)) +
        ylab("normalized counts") +
        ggtitle(paste("mean", strand, "TSS-seq signal"),
                subtitle = paste(nindices, ylabel)) +
        theme_default +
        theme(legend.text = element_text(size=12))
    }
    
    meta_heatmap = function(df){
        ggplot(data = df, aes(x=position, y=0, fill=log2(mean+pcount))) +
        geom_raster() +
        scale_fill_viridis(option="inferno", name=expression(bold(log[2] ~ signal)),
                           guide=guide_colorbar(barheight=8, barwidth=1)) +
        scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                           labels=c(refptlabel, "", endlabel),
                           name="scaled distance",
                           limits = c(-upstream/1000, (dnstream+scaled_length)/1000),
                           expand=c(0,0)) +
        scale_y_continuous(breaks=0,expand=c(0,0), name=NULL) +
        ggtitle(paste("mean", strand, "TSS-seq signal"),
                subtitle = paste(nindices, ylabel)) +
        theme_default +
        theme(axis.text.y = element_blank(),
              strip.text.y = element_text(hjust=1),
              legend.title = element_text(size=10, face="plain"),
              axis.ticks.x = element_line(color="black", size=1),
              axis.ticks.y = element_blank(),
              panel.border = element_blank())
    }
    
    raw = read_tsv(intable, col_names=c("group", "sample", "index", "position","cpm")) %>%
        filter(sample %in% samplelist & !is.na(cpm)) %>% 
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
        
    nindices = max(raw$index, na.rm=TRUE)
    nsamples = length(samplelist)
    ngroups = length(fct_unique(raw$group))

    pcount=0.1
    #get replicate info for sample facetting
    repl_df = raw %>% select(group, sample) %>% distinct() %>% group_by(group) %>%
        mutate(replicate=row_number()) %>% ungroup() %>% select(-group)
    
    #plot heatmap facetted by sample and group
    df_sample = raw %>% group_by(group, sample, position) %>% 
        summarize(mean = winsor.mean(cpm, trim = trim_pct),
                  sem = winsor.sd(cpm, trim = trim_pct)/sqrt(n())) %>% 
        left_join(repl_df, by="sample")
    
    meta_sample = meta(df_sample) +
        facet_grid(replicate~group) +
        theme(strip.text.y = element_text(angle=0))
    ggsave(meta_sample_out, plot=meta_sample, height=2+4.5*max(repl_df$replicate),
           width=7*ngroups, units="cm", limitsize=FALSE)
    rm(meta_sample)
    
    meta_sample_overlay = meta_overlay(df_sample) +
        stat_stepribbon(aes(group=sample)) +
        geom_step(aes(group=sample, color=group),alpha=1)
    ggsave(meta_sample_overlay_out, plot=meta_sample_overlay,
           height=8, width=14, units="cm")
    rm(meta_sample_overlay)
    
    meta_heatmap_sample = meta_heatmap(df_sample) +
        facet_grid(sample~., switch="y")
    ggsave(meta_heatmap_sample_out, plot=meta_heatmap_sample,
           height=2.5+1.25*nsamples, width=14, units="cm")
    rm(meta_heatmap_sample, df_sample, repl_df)
    
    df_group = raw %>% group_by(group, position) %>%
        summarise(mean = winsor.mean(cpm, trim=trim_pct),
                  sem = winsor.sd(cpm, trim=trim_pct)/sqrt(n()))
    rm(raw)
    
    meta_group = meta(df_group) +
        facet_grid(.~group)
    ggsave(meta_group_out, plot=meta_group, height=8,
           width=7*ngroups, units="cm", limitsize=FALSE)
    rm(meta_group_out)
    
    meta_group_overlay = meta_overlay(df_group) +
        stat_stepribbon(aes(group=group)) +
        geom_step(aes(group=group, color=group),alpha=1)
    ggsave(meta_group_overlay_out, plot=meta_group_overlay,
           height=8, width=14, units="cm")
    rm(meta_group_overlay)
    
    meta_heatmap_group = meta_heatmap(df_group) +
        facet_grid(group~., switch="y")
    ggsave(meta_heatmap_group_out, meta_heatmap_group,
           height=2.5+1.5*ngroups, width=14, units="cm")
}

main(intable= snakemake@input[["matrix"]],
     samplelist = snakemake@params[["samplelist"]],
     strand = tolower(snakemake@wildcards[["strand"]]),
     upstream = snakemake@params[["upstream"]],
     dnstream= snakemake@params[["dnstream"]],
     trim_pct = snakemake@params[["trim_pct"]],
     refptlabel = snakemake@params[["refpointlabel"]],
     scaled_length = snakemake@params[["scaled_length"]],
     endlabel = snakemake@params[["endlabel"]],
     ylabel = snakemake@params[["ylabel"]],
     meta_sample_out = snakemake@output[["meta_sample"]],
     meta_sample_overlay_out = snakemake@output[["meta_sample_overlay"]],
     meta_heatmap_sample_out = snakemake@output[["meta_heatmap_sample"]],
     meta_group_out = snakemake@output[["meta_group"]],
     meta_group_overlay_out = snakemake@output[["meta_group_overlay"]],
     meta_heatmap_group_out = snakemake@output[["meta_heatmap_group"]])

