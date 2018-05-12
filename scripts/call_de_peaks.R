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

main = function(intable, norm, sitable, samples, groups, condition, control, alpha, lfc,
                results_all, results_up, results_down, results_unch,
                normcounts, rldcounts, qcplots){
    #import data
    countdata = get_countdata(intable, samples)
    coldata = data.frame(condition=factor(groups,
                                          levels = c(control, condition)),
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
    } else {
        dds = dds %>% estimateSizeFactors()
    }
    dds = dds %>% estimateDispersions() %>% nbinomWaldTest()

    #extract normalized counts and write to file
    ncounts = dds %>%
        counts(normalized=TRUE) %>%
        as.data.frame() %>%
        rownames_to_column(var='name') %>%
        as_tibble()
    ncountsavg = ncounts %>%
        gather(sample, value, -name) %>%
        mutate(group = if_else(sample %in% samples[groups==condition], "condition_expr", "control_expr")) %>%
        group_by(name, group) %>% summarise(mean = mean(value)) %>% spread(group, mean) %>%
        ungroup()

    #plot sd vs. mean for unshrunken (log2) counts
    ntd = dds %>%
        normTransform() %>%
        assay() %>%
        as.data.frame() %>%
        rownames_to_column(var="name") %>%
        as_tibble() %>%
        gather(sample, signal, -name) %>%
        group_by(name) %>%
        summarise(mean=mean(signal), sd=sd(signal)) %>%
        mutate(rank = min_rank(desc(mean)))
    maxsd = max(ntd$sd)*1.01
    ntdplot = mean_sd_plot(ntd, maxsd) +
        ggtitle(expression(paste("raw ", log[2]("counts"))))

    #extract rlog transformed counts and write to file
    rlogcounts = dds %>%
        rlog(blind=FALSE) %>%
        assay() %>%
        as.data.frame() %>%
        rownames_to_column(var="name") %>% as_tibble()
    #plot sd vs. mean for rlog transformed counts
    rld = rlogcounts %>%
        gather(sample, signal, -name) %>%
        group_by(name) %>%
        summarise(mean=mean(signal), sd=sd(signal)) %>%
        mutate(rank = min_rank(desc(mean)))
    rldplot = mean_sd_plot(rld, maxsd) +
        ggtitle(expression(paste("regularized ", log[2]("counts"))))

    #extract DESeq2 results and write to file
    resdf = results(dds, alpha=alpha, lfcThreshold=lfc, altHypothesis="greaterAbs") %>%
        as_data_frame() %>%
        rownames_to_column(var='name') %>%
        arrange(padj) %>%
        inner_join(ncountsavg, by='name') %>%
        rownames_to_column(var="peak_name") %>%
        mutate_at(vars(peak_name), funs(paste0("peak_", .))) %>%
        mutate(score = as.integer(pmin(-125*log2(padj), 1000))) %>%
        mutate_at(c('pvalue','padj'), funs(-log10(.))) %>%
        mutate_if(is.numeric, round, 3) %>%
        dplyr::rename(log2_foldchange=log2FoldChange, lfc_SE=lfcSE,
                      log10_pval=pvalue, log10_padj=padj, mean_expr=baseMean)

    ncounts = resdf %>%
        select(name, peak_name) %>%
        inner_join(ncounts, by='name') %>%
        separate(name, into=c('chrom','strand','start','end'), sep="-") %>%
        mutate_at(vars(strand), funs(if_else(.=="minus", "-", "+"))) %>%
        mutate_if(is.numeric, round, 3) %>%
        mutate_at(vars(start, end), funs(as.integer(.))) %>%
        write_tsv(path=normcounts, col_names=TRUE)

    rlogcounts = resdf %>%
        select(name, peak_name) %>%
        inner_join(rlogcounts, by='name') %>%
        separate(name, into=c('chrom','strand','start','end'), sep="-") %>%
        mutate_at(vars(strand), funs(if_else(.=="minus", "-", "+"))) %>%
        mutate_if(is.numeric, round, 3) %>%
        mutate_at(vars(start, end), funs(as.integer(.))) %>%
        write_tsv(path=rldcounts, col_names=TRUE)

    resdf = resdf %>%
        separate(name, into=c('chrom','strand','start','end'), sep="-") %>%
        mutate_at(vars(strand), funs(if_else(.=="minus", "-", "+"))) %>%
        mutate_if(is.numeric, round, 3) %>%
        mutate_at(vars(start, score, end), funs(as.integer(.))) %>%
        select(chrom, start, end, peak_name, score, strand,
               log2_foldchange, lfc_SE, stat, log10_pval, log10_padj,
               mean_expr, condition_expr, control_expr) %>%
        write_tsv(path=results_all, col_names=TRUE)

    resdf_sig = resdf %>% filter(log10_padj> -log10(alpha))
    resdf_nonsig = resdf %>% filter(log10_padj<= -log10(alpha)) %>%
        write_tsv(results_unch)

    resdf_sig %>%
        filter(log2_foldchange >=0) %>%
        write_tsv(results_up)

    resdf_sig %>%
        filter(log2_foldchange <0) %>%
        write_tsv(results_down)

    maplot = ggplot() +
        geom_hline(yintercept = 0, color="black", linetype="dashed") +
        geom_point(data = resdf_nonsig, aes(x=mean_expr, y=log2_foldchange),
                   color="black", alpha=0.3, stroke=0, size=0.7) +
        geom_point(data = resdf_sig, aes(x=mean_expr, y=log2_foldchange),
                   color="red", alpha=0.3, stroke=0, size=0.7) +
        scale_x_log10(name="mean of normalized counts") +
        ylab(substitute(log[2]~frac(cond,cont), list(cond=condition, cont=control))) +
        theme_light() +
        theme(text = element_text(size=8))

    volcano = ggplot() +
        geom_point(data = resdf_nonsig, aes(x=log2_foldchange, y = log10_padj),
                   alpha=0.3, stroke=0, size=0.7) +
        geom_point(data = resdf_sig, aes(x=log2_foldchange, y = log10_padj),
                   alpha=0.3, stroke=0, size=0.7) +
        geom_hline(yintercept = -log10(alpha), color="red", linetype="dashed") +
        xlab(substitute(log[2]~frac(cond,cont), list(cond=condition, cont=control))) +
        ylab(expression(-log[10]("p value"))) +
        theme_light() +
        theme(text = element_text(size=8))

    #plot library size vs sizefactor
    sfdf = dds %>%
        sizeFactors() %>%
        as_tibble() %>%
        rownames_to_column(var="sample") %>%
        dplyr::rename(sizefactor=value) %>%
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

    out = arrangeGrob(sfplot, ggplot()+theme_void(),
                      ntdplot, rldplot,
                      maplot, volcano, ncol=2,
                      heights = unit(c(4, 6, 6), rep("cm",3)))
    ggsave(qcplots, out, height=18, width = 16, units="cm")
}

main(intable = snakemake@input[["expcounts"]],
     norm = snakemake@wildcards[["norm"]],
     sitable = snakemake@input[["sicounts"]],
     samples = snakemake@params[["samples"]],
     groups = snakemake@params[["groups"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     alpha = snakemake@params[["alpha"]],
     lfc = snakemake@params[["lfc"]],
     results_all = snakemake@output[["results_all"]],
     results_up = snakemake@output[["results_up"]],
     results_down = snakemake@output[["results_down"]],
     results_unch = snakemake@output[["results_unch"]],
     normcounts = snakemake@output[["normcounts"]],
     rldcounts = snakemake@output[["rldcounts"]],
     qcplots = snakemake@output[["qcplots"]])

