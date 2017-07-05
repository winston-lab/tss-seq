library(tidyverse)
library(forcats)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(viridis)
library(GGally)

import = function(path){
  read_table2(path, col_names=TRUE) 
}

get_countdata = function(table, samplenames){
  df = data.frame(table[,-1], row.names=table$'chrom-start-end')
  names(df) = samplenames
  return(df)
}

extract_deseq_results = function(dds, alpha){
  results(dds, alpha=alpha) %>% as.data.frame() %>% rownames_to_column(var='base') %>% as_data_frame() %>% filter(padj != "NA") %>% filter(padj < alpha) %>% arrange(padj) %>% return()
}

plot_correlation = function(path, dds){
  df = dds %>% counts(normalized=TRUE) %>% as.data.frame()
  lim = df %>% log10() %>% max() %>% ceiling()

  ggsave(path, plot=
    ggscatmat(data = df, alpha = 0, corMethod = "pearson") +
              geom_point(size=.001, alpha=0.4) +
              scale_x_log10(limits = c(1, 1*10^lim)) +
              scale_y_log10(limits = c(1, 1*10^lim)) +
              coord_fixed(ratio=1) +
              theme_bw() +
              xlab("normalized counts") +
              ylab("normalized counts") +
              theme(axis.text = element_text(size=12),
                    axis.title = element_text(size=12),
                    strip.text = element_text(size=12, face="bold"),
                    strip.background = element_blank()),
    width = 5*ncol(df),
    height = 4*ncol(df),
    units = "cm")
}

plot_count_heatmap = function(path, counts, de_results, alpha){
  data = assay(counts)[row.names(assay(counts)) %in% (de_results %>% filter(padj < alpha))$base,]
  pheatmap(data,
           cluster_rows = TRUE,
           show_rownames = FALSE,
           cluster_cols = FALSE,
           col = viridis(100),
           border_color=NA,
           width = ncol(assay(counts)),
           height = 4+.01*nrow(data),
           fontsize=12,
           filename=path,
           main = paste("rlog-transformed counts for\n", nrow(data), "DE bases @ FDR =", alpha))
}

plot_dist_heatmap = function(path, counts){
  sampleDists = counts %>% assay() %>% t() %>% dist()
  sampleDistMatrix = sampleDists %>% as.matrix()
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = inferno(100),
           border_color = NA,
           fontsize = 12,
           main = "sample-to-sample Euclidean distances of\n rlog-transformed counts",
           cellheight = 50,
           cellwidth = 50,
           filename = path)
}

plot_pca = function(pcapath, screepath, counts){
  pca = counts %>% assay() %>% t() %>% prcomp(center=TRUE, scale=FALSE)
  pca.df = data.frame(pca$x, group = colData(counts)[["condition"]]) %>% rownames_to_column(var="name")
  percentVar = data.frame(pctvar = (pca$sdev^2)/sum(pca$sdev^2)) %>% rownames_to_column(var="PC") %>% as_tibble()
  percentVar$PC = percentVar$PC %>% fct_inorder()
  
  scree = ggplot(data = percentVar, aes(x=PC, y=pctvar*100)) +
                  geom_col() +
                  theme_minimal() +
                  ylab("% variance explained") +
                  theme(axis.title = element_text(size=12, face="bold"),
                        axis.text = element_text(size=12))
  
  p = ggplot(data = pca.df, aes(x=PC1, y=PC2)) +
              geom_point(aes(color=group), size=4) +
              geom_text_repel(aes(label=name), size=3, point.padding = unit(.2, "cm")) +
              scale_color_brewer(palette = "Set1", direction=-1) +
              xlab(paste("PC1: ", round(percentVar[1,2]*100, 1), "% of variance explained", sep = "")) +
              ylab(paste("PC2: ", round(percentVar[2,2]*100, 1), "% \nof variance explained", sep = "")) +
              theme_minimal() +
              theme(legend.position="none",
                    axis.text = element_text(size=12),
                    axis.title = element_text(size=12, face="bold"))
  
  ggsave(pcapath, p, width = 12, height = max(12*percentVar[2,2]/percentVar[1,2], 6), units = "cm")
  ggsave(screepath, scree, width = nrow(percentVar), height = 6, units = "cm")
}

plot_size_v_sf = function(path, title, dds){
  sf = dds$sizeFactor %>% as.data.frame() %>% rownames_to_column(var="sample") %>% as_tibble()
  names(sf) = c("sample", "sizefactor")
  df = dds %>% counts() %>% as.data.frame() %>% rownames_to_column(var="base") %>%
    gather(key=sample, value=counts, -base) %>% group_by(sample) %>% summarise(totalcounts = sum(counts)) %>%
    left_join(sf, by="sample") 
  
  p = ggplot(data = df, aes(x=totalcounts/1e6, y=sizefactor)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE, color="red") +
    geom_text_repel(aes(label = sample), size=4) +
    xlab("library size (million reads)") +
    ylab("DESeq2 size factor") +
    ggtitle(title) +
    theme_minimal() +
    theme(axis.title = element_text(size=12, face="bold"),
          axis.text = element_text(size=12),
          title = element_text(size=12, face="bold"))
  
  ggsave(path, plot = p, width = 11, height = 7, units = "cm")
  return(df)
}

plot_spikein_pct = function(sfpath1, sfpath2, sipctpath, dds, si.dds){
  df = plot_size_v_sf(path = sfpath1, title = "experimental", dds)
  df.si = plot_size_v_sf(path = sfpath2, title = "spike-in", si.dds)
  
  names(df.si) = c("sample", "spikecounts", "spikesizefactor")
  
  df = df %>% left_join(df.si, by="sample") %>% mutate(si.pct = spikecounts/(spikecounts+totalcounts))
  
  p = ggplot(data = df, aes(x=sample, y = round(si.pct*100, 2))) +
    geom_col() +
    geom_text(aes(label = round(si.pct*100, 1)), size=5, position = position_stack(vjust = .9)) +
    theme_minimal() +
    ylab("% spike-in\nreads") +
    theme(axis.title = element_text(size=12, face="bold"),
          axis.title.y = element_text(angle=0, vjust = 0.5),
          axis.title.x = element_blank(),
          axis.text = element_text(size=12, color="black"),
          axis.text.x = element_text(angle=30, hjust = 0.8))
  
  ggsave(sipctpath, plot = p, width = 1.5+2*nrow(df), height = 6, units = "cm")
  
  #ggplot(data = df, aes(x=sample, y = spikesizefactor)) +
  #        geom_col() +
  #        geom_text(aes(label = round(spikesizefactor, 3)), size=5, position = position_stack(vjust = .9)) +
  #        theme_minimal()
  
}

qual_ctrl = function(intable,
                     intable.si,
                     samplenames,
                     samplegroups,
                     nospikein,
                     corrplot.spikenorm,
                     corrplot.libsizenorm,
                     alpha,
                     #count.heatmap.spikenorm,
                     #count.heatmap.libsizenorm,
                     dist.heatmap.spikenorm,
                     dist.heatmap.libsizenorm,
                     pca.spikenorm,
                     scree.spikenorm,
                     pca.libsizenorm,
                     scree.libsizenorm,
                     de.spikenorm.path,
                     de.libsizenorm.path){
  raw = import(intable)
  countdata = get_countdata(raw, samplenames = samplenames)
  if(is.null(nospikein)){
      raw.si = import(intable.si)
      countdata.si = get_countdata(raw.si, samplenames = samplenames) 
  }
  coldata.all = data.frame(condition=factor(samplegroups, levels = unique(samplegroups)), row.names=names(countdata))
  
  dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata.all, design = ~condition)
  if(is.null(nospikein)){
      si.dds = DESeqDataSetFromMatrix(countData = countdata.si, colData = coldata.all, design= ~condition)
      si.dds = estimateSizeFactors(si.dds)
      dds.spikenorm = dds 
      sizeFactors(dds.spikenorm) = sizeFactors(si.dds)
      dds.spikenorm = dds.spikenorm %>% estimateDispersions() %>% nbinomWaldTest()
      
      plot_correlation(corrplot.spikenorm, dds.spikenorm)
      resdf = extract_deseq_results(dds.spikenorm, alpha=alpha)
      write.table(resdf, file=de.spikenorm.path, quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
      rld = rlog(dds.spikenorm, blind=FALSE)
      rld.df = rld %>% assay() %>% as.data.frame() %>% rownames_to_column() %>% as_data_frame()
  
      #plot_count_heatmap(count.heatmap.spikenorm, rld, resdf, alpha)
      plot_dist_heatmap(dist.heatmap.spikenorm, rld)
      plot_pca(pca.spikenorm, scree.spikenorm, rld)
  } else{
         file.create(corrplot.spikenorm,
                        #count.heatmap.spikenorm,
                        dist.heatmap.spikenorm,
                        pca.spikenorm,
                        scree.spikenorm,
                        de.spikenorm.path)
  }
  dds  = dds %>% estimateSizeFactors() %>% estimateDispersions() %>% nbinomWaldTest()
  
  plot_correlation(corrplot.libsizenorm, dds)
  resdf.nospike = extract_deseq_results(dds, alpha=alpha)
  write.table(resdf.nospike, file=de.libsizenorm.path, quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
  rld.nospike = rlog(dds, blind=FALSE)
  rld.nospike.df = rld.nospike %>% assay() %>% as.data.frame() %>% rownames_to_column() %>% as_data_frame()
  #plot_count_heatmap(count.heatmap.libsizenorm, rld.nospike, resdf.nospike, alpha)
  plot_dist_heatmap(dist.heatmap.libsizenorm, rld.nospike)
  plot_pca(pca.libsizenorm, scree.libsizenorm, rld.nospike)
}

qc = qual_ctrl(intable = snakemake@input[["exp"]],
                    intable.si = snakemake@input[["si"]],
                    samplenames = snakemake@params[["samples"]],
                    samplegroups = snakemake@params[["samplegroups"]],
                    nospikein = snakemake@params[["nospikein"]],
                    corrplot.spikenorm = snakemake@output[["corrplot_spikenorm"]],
                    corrplot.libsizenorm = snakemake@output[["corrplot_libsizenorm"]],
                    alpha=snakemake@params[["alpha"]],
                    #count.heatmap.spikenorm = snakemake@output[["count_heatmap_spikenorm"]],
                    #count.heatmap.libsizenorm = snakemake@output[["count_heatmap_libsizenorm"]],
                    dist.heatmap.spikenorm = snakemake@output[["dist_heatmap_spikenorm"]],
                    dist.heatmap.libsizenorm = snakemake@output[["dist_heatmap_libsizenorm"]],
                    pca.spikenorm = snakemake@output[["pca_spikenorm"]],
                    scree.spikenorm = snakemake@output[["scree_spikenorm"]],
                    pca.libsizenorm = snakemake@output[["pca_libsizenorm"]],
                    scree.libsizenorm = snakemake@output[["scree_libsizenorm"]],
                    de.spikenorm.path = snakemake@output[["de_spikenorm_path"]],
                    de.libsizenorm.path = snakemake@output[["de_libsizenorm_path"]])
