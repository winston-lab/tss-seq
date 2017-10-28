library(tidyverse)
library(forcats)
library(DESeq2)
library(viridis)
library(scales)
library(gridExtra)
library(ggrepel)

get_countdata = function(path, samples){
    df = read_tsv(path, col_names=TRUE) %>% select(c("name", samples)) %>% 
            column_to_rownames(var="name") %>% as.data.frame()
    df = df[rowSums(df)>1,]
    return(df)
}

extract_deseq_results = function(dds, alpha){
    results(dds, alpha=alpha) %>% as_data_frame() %>%
        rownames_to_column(var="name") %>% arrange(padj) %>% return()
}

mean_sd_plot = function(df, ymax){
    ggplot(data = df, aes(x=rank, y=sd)) +
        geom_hex(aes(fill=log10(..count..), color=log10(..count..)), bins=40, size=0) +
        geom_smooth(color="#4292c6") +
        scale_fill_viridis(option="inferno", name=expression(log[10](count)), guide=FALSE) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        scale_x_continuous(trans=reverselog_trans(10), name="rank(mean expression)") +
        scale_y_continuous(limits = c(NA, ymax)) +
        ylab("SD") +
        theme_light() +
        theme(text = element_text(size=8))
}

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

call_de_bases = function(intable, norm, sitable, samples, groups, condition, control, alpha, results, normcounts, rldcounts, qcplots){
    #import data 
    countdata = get_countdata(intable, samples)
    coldata = data.frame(condition=factor(groups,
                                          levels = unique(groups)),
                         row.names=names(countdata))
    #run DESeq2 
    dds = DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ condition)
    
    if (norm=="spikenorm"){
        sicountdata = get_countdata(sitable, samples)
        sidds = DESeqDataSetFromMatrix(countData = sicountdata,
                                       colData = coldata,
                                       design = ~ condition)
        sidds = sidds %>% estimateSizeFactors()
        sizeFactors(dds) = sizeFactors(sidds)
    }
    else {
        dds = dds %>% estimateSizeFactors()
    }
    dds = dds %>% estimateDispersions() %>% nbinomWaldTest()
    
    #extract DESeq2 results and write to file
    resdf = results(dds, alpha=alpha) %>% as_data_frame() %>%
                rownames_to_column(var="name") %>% arrange(padj) 
    write_tsv(resdf, path=results, col_names=TRUE)
    
    #extract normalized counts and write to file
    ncounts = dds %>% counts(normalized=TRUE) %>% as.data.frame() %>%
                rownames_to_column(var="name") %>% as_tibble()
    write_tsv(ncounts, path=normcounts, col_names=TRUE)
    
    #plot sd vs. mean for unshrunken (log2) counts
    ntd = dds %>% normTransform() %>% assay() %>% as.data.frame() %>%
            rownames_to_column(var="name") %>% as_tibble() %>%
            gather(sample, signal, -name) %>% group_by(name) %>%
            summarise(mean=mean(signal), sd=sd(signal)) %>%
            mutate(rank = min_rank(desc(mean)))
    maxsd = max(ntd$sd)*1.01
    ntdplot = mean_sd_plot(ntd, maxsd) + 
                ggtitle(expression(paste("raw ", log[2]("counts"))))
        
    #extract rlog transformed counts and write to file
    rld = dds %>% rlog(blind=FALSE) %>% assay() %>% as.data.frame() %>%
            rownames_to_column(var="name") %>% as_tibble()
    write_tsv(rld, path=rldcounts, col_names=TRUE)
    
    #plot sd vs. mean for rlog transformed counts
    rld = rld %>% gather(sample, signal, -name) %>% group_by(name) %>%
            summarise(mean=mean(signal), sd=sd(signal)) %>%
            mutate(rank = min_rank(desc(mean)))
    rldplot = mean_sd_plot(rld, maxsd) +
                ggtitle(expression(paste("regularized ", log[2]("counts"))))
    
    #plot library size vs sizefactor
    sfdf = dds %>% sizeFactors() %>% as_tibble() %>%
            rownames_to_column(var="sample") %>% dplyr::rename(sizefactor=value) %>%
            inner_join(colSums(countdata) %>% as_tibble() %>%
                       rownames_to_column(var="sample") %>% dplyr::rename(libsize=value),
                       by="sample")
    sfplot = ggplot(data = sfdf, aes(x=libsize/1e6, y=sizefactor)) +
                geom_smooth(method="lm", se=FALSE, color="red", size=0.5) +
                geom_point(size=0.5) +
                geom_text_repel(aes(label=sample), size=2) +
                xlab("library size (M reads)") +
                ylab("size factor (median of ratios)") +
                theme_light() +
                theme(text = element_text(size=8))
    
    #MA plot for differential expression
    resdf.sig = resdf %>% filter(padj<alpha)
    resdf.nonsig = resdf %>% filter(padj>=alpha)
    maplot = ggplot() +
                geom_hline(yintercept = 0, color="black", linetype="dashed") +
                geom_point(data = resdf.nonsig, aes(x=baseMean, y=log2FoldChange),
                           color="black", alpha=0.3, stroke=0, size=0.7) +
                geom_point(data = resdf.sig, aes(x=baseMean, y=log2FoldChange),
                           color="red", alpha=0.3, stroke=0, size=0.7) +
                scale_x_log10(name="mean of normalized counts") +
                ylab(substitute(log[2]~frac(cond,cont), list(cond=condition, cont=control))) +
                theme_light() +
                theme(text = element_text(size=8))
    
    volcano = ggplot() +
                geom_point(data = resdf.nonsig, aes(x=log2FoldChange, y = -log10(padj)),
                           alpha=0.3, stroke=0, size=0.7) +
                geom_point(data = resdf.sig, aes(x=log2FoldChange, y = -log10(padj)),
                           alpha=0.3, stroke=0, size=0.7) +
                geom_hline(yintercept = -log10(alpha), color="red", linetype="dashed") +
                xlab(substitute(log[2]~frac(cond,cont), list(cond=condition, cont=control))) +
                ylab(expression(-log[10]("p value"))) +
                theme_light() +
                theme(text = element_text(size=8))
    
    out = grid.arrange(sfplot, ggplot()+theme_void(),
                       ntdplot, rldplot,
                       maplot, volcano, ncol=2,
                       heights = unit(c(4, 6, 6), rep("cm",3)))
    ggsave(qcplots, out, height=18, width = 16, units="cm")
}

qc = call_de_bases(intable = snakemake@input[["expcounts"]],
                   norm = snakemake@wildcards[["norm"]],
                   sitable = snakemake@input[["sicounts"]],
                   samples = snakemake@params[["samples"]],
                   groups = snakemake@params[["groups"]],
                   condition = snakemake@wildcards[["condition"]],
                   control = snakemake@wildcards[["control"]],
                   alpha = snakemake@params[["alpha"]],
                   results = snakemake@output[["results"]],
                   normcounts = snakemake@output[["normcounts"]],
                   rldcounts = snakemake@output[["rldcounts"]],
                   qcplots = snakemake@output[["qcplots"]])
