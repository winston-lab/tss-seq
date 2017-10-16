library(tidyverse)
library(forcats)
library(viridis)

plotheatmaps = function(intable, upstream, dnstream, pct_cutoff, refptlab, ylabel, cmap, samples_out, group_out){
    raw = read_tsv(intable,
    	 col_names=c("group", "sample", "index", "position","cpm"),
    	 col_types=cols(group=col_character(), sample=col_character(), index=col_integer(), position=col_double(), cpm=col_double())) %>% filter(cpm != "NA")
    raw$sample = fct_inorder(raw$sample, ordered = TRUE)
    raw$group = fct_inorder(raw$group, ordered = TRUE)
    
    nindices = max(raw$index, na.rm=TRUE)
    nsamples = length(fct_unique(raw$sample))
    ngroups = length(fct_unique(raw$group))
    #w = round((max(raw$position) - min(raw$position))*1000/148)
    
    #pseudocount for log-transform
    pcount = .01
    
    #percentile cutoff for heatmap visualization
    cutoff = quantile(raw$cpm, probs=pct_cutoff, na.rm=TRUE)
    
    #plot heatmap facetted by sample and group
    heatmap_base = ggplot(data = raw %>% mutate_at(vars(cpm), funs(pmin(cutoff, .)))) +
      geom_raster(aes(x=position, y=index, fill=log2(cpm+pcount))) +
      scale_y_reverse(name=paste(nindices, ylabel), expand=c(0.01, 0))+
      scale_x_continuous(breaks = c(-upstream/1000, 0, dnstream/1000),
                         labels=c(ifelse(upstream>200, -upstream/1000, ''),
                                  refptlab,
                                  ifelse(dnstream>200, dnstream/1000, '')),
                         minor_breaks = scales::pretty_breaks(n=10),
                         name=paste("distance from", refptlab, "(kb)")) +
      scale_fill_viridis(option = cmap, na.value="#FFFFFF00", name=expression(paste(log[2], '(TSS-seq signal)')), guide=guide_colorbar(title.position="top", barwidth=15, barheight=1, title.hjust=0.5)) +
      theme_minimal() +
      theme(text = element_text(size=12, face="bold", color="black"),
              legend.position = "top",
              legend.text = element_text(size=8, face="plain"),
              strip.text = element_text(size=12, face="bold", color="black"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=12, face="bold", color="black", margin = unit(c(0,0,0,0),"cm")),
              panel.grid.major.x = element_line(color="black"),
              panel.grid.minor.x = element_line(color="grey80"),
              panel.grid.major.y = element_line(color="grey80"),
              panel.grid.minor.y = element_blank(),
              panel.spacing.x = unit(.5, "cm"))
    
    hmap.width = max(12, ((upstream+dnstream)/200)*ngroups)

    heatmap_samples = heatmap_base + facet_wrap(~sample, dir="v", ncol=ngroups)
    ggsave(samples_out , plot = heatmap_samples,
           height= (.0005*nindices+7.5)*nsamples/ngroups,
           width = hmap.width, units = "cm")
    rm(heatmap_samples)
    gc()
    heatmap_groups = heatmap_base + facet_wrap(~group, ncol=ngroups)
    ggsave(group_out, plot = heatmap_groups,
           height= .002*nindices+14.75,
           width = hmap.width, units = "cm")
}

plotheatmaps(intable= snakemake@input[["matrix"]],
             upstream = snakemake@params[["upstream"]],
             dnstream= snakemake@params[["dnstream"]],
             pct_cutoff = snakemake@params[["pct_cutoff"]],
             refptlab = snakemake@params[["refpointlabel"]],
             ylabel = snakemake@params[["ylabel"]],
             cmap = snakemake@params[["heatmap_cmap"]],
             samples_out = snakemake@output[["heatmap_sample"]],
             group_out = snakemake@output[["heatmap_group"]])
