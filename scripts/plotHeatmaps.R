library(tidyverse)
library(forcats)
library(viridis)

raw = read_table2(snakemake@input[["matrix"]],
	 col_names=c("group", "sample", "index", "position","cpm"),
	 col_types=cols(group=col_character(), sample=col_character(), index=col_integer(), position=col_double(), cpm=col_double())) %>% filter(cpm != "NA")
raw$sample = factor(raw$sample, ordered = TRUE)
raw$group = factor(raw$group, ordered = TRUE)

nindices =  max(raw$index)
nsamples = length(fct_unique(raw$sample))
ngroups = length(fct_unique(raw$group))
w = round((max(raw$position) - min(raw$position))*1000/148)

upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["dnstream"]]

#pseudocount for log-transform
pcount = .01

#percentile cutoff for heatmap visualization
#cutoff = quantile(raw$cpm, probs=snakemake@params[["pct_cutoff"]], na.rm=TRUE)

#plot heatmap facetted by sample and group
heatmap_base = ggplot(data = raw) +
  geom_tile(aes(x=position, y=index, fill=log2(cpm+pcount))) +
  scale_y_reverse(name=paste(nindices, snakemake@params[["ylabel"]])) +
  scale_x_continuous(breaks = c(-upstream/1000, 0, downstream/1000), labels=c(-upstream/1000, snakemake@params[["refpointlabel"]], downstream/1000)) +
  xlab(paste("distance from", snakemake@params[["refpointlabel"]], "(kb)")) +
  scale_fill_viridis(option = snakemake@params[["heatmap_cmap"]],na.value="white", name=expression(paste(log[2], '(TSS-seq normalized counts)')), guide=guide_colorbar(title.position="top", barwidth=15, barheight=1, title.hjust=0.5)) +
  theme_minimal() +
  theme(strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.ticks.length = unit(-2, "mm"))

heatmap_samples = heatmap_base + facet_wrap(~sample, ncol=(nsamples/ngroups))
ggsave(snakemake@output[["heatmap_sample"]], plot = heatmap_samples, height=3+round((nindices/1000)*(ngroups)), width = 3+.3*w*(nsamples/ngroups), units = "cm")
rm(heatmap_samples)
heatmap_groups = heatmap_base + facet_wrap(~group)
ggsave(snakemake@output[["heatmap_group"]], plot = heatmap_groups, height=3+round(nindices/600), width = 3+.3*w*ngroups, units = "cm")
