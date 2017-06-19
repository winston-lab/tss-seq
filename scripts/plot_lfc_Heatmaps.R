library(tidyverse)
library(forcats)
library(pals)
library(scales)

raw = read_table2(snakemake@input[["matrix"]],
	 col_names=c("controlgroup", "conditiongroup", "control", "condition", "index", "position","lfc"),
	 col_types=cols(controlgroup=col_character(), conditiongroup=col_character(), control=col_character(), condition=col_character(), index=col_integer(), position=col_double(), lfc=col_double())) %>% filter(lfc != "NA")
raw$control = factor(raw$control, ordered = TRUE)
raw$condition = factor(raw$condition, ordered = TRUE)
raw$controlgroup= factor(raw$controlgroup, ordered = TRUE)
raw$conditiongroup= factor(raw$conditiongroup, ordered = TRUE)

nindices =  max(raw$index)
nconditions = length(fct_unique(raw$condition))
ncontrols= length(fct_unique(raw$control))
nconditiongroups = length(fct_unique(raw$conditiongroup))
ncontrolgroups= length(fct_unique(raw$controlgroup))

w = round((max(raw$position) - min(raw$position))*1000/148)

upstream = snakemake@params[["upstream"]]
downstream = snakemake@params[["dnstream"]]

#plot heatmap facetted by sample and group
heatmap_base = ggplot(data = raw) +
  geom_tile(aes(x=position, y=index, fill=lfc)) +
  scale_y_reverse(name=paste(nindices, snakemake@params[["ylabel"]])) +
  scale_x_continuous(breaks = c(-upstream/1000, 0, downstream/1000), labels=c(-upstream/1000, snakemake@params[["refpointlabel"]], downstream/1000)) +
  xlab(paste("distance from", snakemake@params[["refpointlabel"]], "(kb)")) +
  scale_fill_gradientn(colors = kovesi.diverging_bky_60_10_c30(200), values = rescale(c(min(raw$lfc), 0, max(raw$lfc))) , guide=guide_colorbar(title.position="top", barwidth = 15, barheight=1, title.hjust = 0.5), name=expression(paste(log[2], bgroup("(", frac("condition TSS-seq counts", "control TSS-seq counts"), ")")))) +
  theme_minimal() +
  theme(strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.ticks.length = unit(-2, "mm"))

heatmap_samples = heatmap_base + facet_grid(condition~control)
ggsave(snakemake@output[["heatmap_sample"]], plot = heatmap_samples, height=3+round((nindices/600)*nconditions), width = 3+.3*w*ncontrols, units = "cm")
rm(heatmap_samples)
heatmap_groups = heatmap_base + facet_grid(conditiongroup~controlgroup)
ggsave(snakemake@output[["heatmap_group"]], plot = heatmap_groups, height=3+round(3*(nindices/1000))*nconditiongroups, width = 3+.3*w*ncontrolgroups, units = "cm")
