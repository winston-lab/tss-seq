library(tidyverse)
library(GGally)
library(viridis)
library(pheatmap)
library(ggrepel)

raw = read_table2(snakemake@input[[1]], guess_max=4000)
window = as.integer(raw[1,3]-raw[1,2])
cpm = raw[,-c(1:3)]
cpm = cpm[which(rowSums(cpm)>0), order(names(cpm))]
n = ncol(cpm)

#pairwise scatter plots
ggsave(snakemake@output[["scatter"]], plot= ggscatmat(data = as.data.frame(cpm), alpha = 0.025, corMethod = "pearson") +
                          scale_x_log10() +
                          scale_y_log10() +
                          coord_fixed(ratio=1) +
                          theme_bw() +
                          xlab(paste("counts per million reads, in", window, "bp non-overlapping windows")) +
                          ylab(paste("CPM,", window, "bp windows")) +
                          theme(axis.text = element_text(size=20),
                                axis.title = element_text(size=20),
                                strip.text = element_text(size=20, face="bold"),
                                strip.background = element_blank()),
         height=2+8*n, width=2+8*n, units="cm")

#heatmap of sample-to-sample distances, clustered and unclustered
sampledists = dist(t(cpm))
sampledistmat = as.matrix(sampledists)
pheatmap(sampledistmat, 
         col=inferno(255),
         main=paste("sample-to-sample\nEuclidean distances"),
         border_color=NA,
         cellwidth = 40,
         cellheight = 40,
         clustering_distance_rows = sampledists,
         clustering_distance_cols = sampledists,
         fontsize = 12,
         width = 1.5+1.5*n,
         height = 1+1.5*n,
         filename = snakemake@output[["dists_cluster"]])


pheatmap(sampledistmat, 
         col=inferno(255),
         main=paste("sample-to-sample\nEuclidean distances"),
         border_color=NA,
         cellwidth = 40,
         cellheight = 40,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 12,
         width = 1.5+1.5*n,
         height = 1.5*n,
         filename = snakemake@output[["dists_nocluster"]])

#pca, plot PC1 vs PC2 and scree plot
cpm.pca = prcomp(t(cpm), center=TRUE, scale=FALSE)
cpm.pca.df = as.data.frame(cpm.pca$x) %>% rownames_to_column(var="sample") %>% as_data_frame()
eig = (cpm.pca$sdev)^2
varpct = eig*100/sum(eig)
cumvar = cumsum(varpct)

ggsave(snakemake@output[["pca"]], plot=ggplot(data = cpm.pca.df, aes(x=PC1, y=PC2)) +
                        geom_point() +
                        geom_text_repel(aes(label = sample), size=5) +
                        xlab(paste("PC1: ", round(varpct[1],1), "% of variance explained", sep="")) +
                        ylab(paste("PC2: ", round(varpct[2],1), "% of variance explained", sep="")) +
                        theme_bw() +
                        theme(axis.text = element_text(size=12)),
       height=10, width=10, units = "cm")


varpct.df = as.data.frame(varpct) %>% rownames_to_column(var="PC")
npc = nrow(varpct.df)

ggsave(snakemake@output[["scree"]], plot = ggplot(data = varpct.df, aes(x=PC, y=varpct)) +
                            geom_col() +
                            xlab("PC") +
                            ylab("% of variance explained") +
                            theme_bw() +
                            theme(axis.text = element_text(size=12),
                                  axis.title = element_text(size=12)),
       height=6, width = 3+npc, units = "cm")
