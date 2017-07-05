library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(viridis)
library(GGally)

get_countdata = function(table, samplenames){
    df = data.frame(table[,-1], row.names=table$X1)
    names(df) = samplenames
    return(df)
}

extract_deseq_results = function(dds, alpha){
    results(dds, alpha=alpha) %>% as.data.frame() %>% rownames_to_column(var='base') %>% as_data_frame() %>% filter(padj != "NA") %>% filter(padj<alpha) %>% arrange(padj) %>% return()
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
               xlab("normalized counts per TSS cluster") +
               ylab("normalized counts per TSS cluster") +
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
             height = max(8, .0003*nrow(data)),
             fontsize=12,
             filename=path,
             main = paste("rlog-transformed counts for\n", nrow(data), "DE TSS clusters @ FDR =", alpha))
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
             main = "sample-to-sample Euclidean distances of\nrlog-transformed counts over TSS clusters",
             cellheight = 50,
             cellwidth = 50,
             filename = path)
}

plot_pca = function(pcapath, screepath, counts){
    pca = counts %>% assay() %>% t() %>% prcomp(center=TRUE, scale=FALSE)
    pca.df = data.frame(pca$x, group = colData(counts)[["condition"]]) %>% rownames_to_column(var="name")
    percentVar = data.frame(pctvar = (pca$sdev^2)/sum(pca$sdev^2)) %>% rownames_to_column(var="PC")
    
    scree = ggplot(data = percentVar, aes(x=PC, y=pctvar*100)) +
        geom_col() +
        theme_minimal() +
        ylab("% variance explained") +
        theme(axis.title = element_text(size=12, face="bold"),
              axis.text = element_text(size=12))
    
    p = ggplot(data = pca.df, aes(x=PC1, y=PC2)) +
        geom_point(aes(color=group), size=4) +
        geom_text_repel(aes(label=name), size=4, point.padding = unit(.2, "cm")) +
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

qual_ctrl = function(intable.cluster,
                     intable.lib,
                     samplenames,
                     samplegroups,
                     corrplot,
                     alpha,
                     count.heatmap,
                     dist.heatmap,
                     pca,
                     scree,
                     de.path){
    
    raw.clusters = read_table2(intable.cluster, col_names=FALSE)
    raw.lib = read_table2(intable.lib, col_names=TRUE)
    names(raw.lib)[1] = "X1" 
    countdata.spikenorm = get_countdata(raw.clusters, samplenames = samplenames)
    #countdata.libsizenorm = get_countdata(raw.libsizenorm, samplenames = samplenames)
    countdata.si = get_countdata(raw.lib, samplenames = samplenames)
    
    coldata = data.frame(condition=factor(samplegroups, levels = unique(samplegroups)), row.names=names(countdata.spikenorm))
    
    dds.clusters = DESeqDataSetFromMatrix(countData = countdata.spikenorm, colData = coldata, design = ~condition)
    dds.lib = DESeqDataSetFromMatrix(countData = countdata.si, colData = coldata, design= ~condition)
    
    #get size factors from spike-in
    dds.lib = estimateSizeFactors(dds.lib)
    sizeFactors(dds.clusters) = sizeFactors(dds.lib)
    
    #do differential expression +/- spike-in
    dds.clusters = dds.clusters %>% estimateDispersions() %>% nbinomWaldTest()
    
    plot_correlation(corrplot, dds.clusters)
    
    resdf = extract_deseq_results(dds.clusters, alpha=alpha)
    write.table(resdf, file=de.path, quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
    
    #transformations for datavis and quality control
    #blinding is FALSE for datavis purposes
    rld = rlog(dds.clusters, blind=FALSE)
    rld.df = rld %>% assay() %>% as.data.frame() %>% rownames_to_column() %>% as_data_frame()
    
    plot_count_heatmap(count.heatmap, rld, resdf, alpha)
    plot_dist_heatmap(dist.heatmap, rld)
    plot_pca(pca, scree, rld)
}

qc = qual_ctrl(intable.cluster = snakemake@input[["clustercounts"]],
               intable.lib = snakemake@input[["libcounts"]],
               samplenames = snakemake@params[["samples"]],
               samplegroups = snakemake@params[["samplegroups"]],
               corrplot = snakemake@output[["corrplot"]],
               alpha= snakemake@params[["alpha"]],
               count.heatmap = snakemake@output[["count_heatmap"]],
               dist.heatmap = snakemake@output[["dist_heatmap"]],
               pca = snakemake@output[["pca"]],
               scree = snakemake@output[["scree"]],
               de.path = snakemake@output[["de_path"]])

