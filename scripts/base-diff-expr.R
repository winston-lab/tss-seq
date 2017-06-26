library(tidyverse)
library(DESeq2)

get_de_bases = function(alpha, samplenames, samplegroups, intable, si.intable, outtable.spikenorm, outtable.libsizenorm){

    raw = read_table2(intable, col_names=FALSE)
    raw.si = read_table2(si.intable, col_names=FALSE)

    countdata = data.frame(raw[,-1], row.names=raw$X1)
    countdata.si = data.frame(raw.si[,-1], row.names=raw.si$X1)
    
    names(countdata) = names(countdata.si) = samplenames

    coldata = data.frame(condition=factor(samplegroups, levels = unique(samplegroups)), row.names=names(countdata))
    
    dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
    si.dds = DESeqDataSetFromMatrix(countData = countdata.si, colData = coldata, design= ~condition)

    #get size factors from spike-in
    si.dds = estimateSizeFactors(si.dds)
    dds.nospike = dds
    sizeFactors(dds) = sizeFactors(si.dds)
 
    #complete DESeq standard workflow +/- spike-in
    dds = estimateDispersions(dds)
    dds = nbinomWaldTest(dds)

    dds.nospike = estimateSizeFactors(dds.nospike)
    dds.nospike = estimateDispersions(dds.nospike)
    dds.nospike = nbinomWaldTest(dds.nospike)

    #get results
    resdf = results(dds, alpha = alpha) %>% as.data.frame() %>% rownames_to_column(var='base') %>% as_data_frame()
    resdf.nospike = results(dds.nospike, alpha = alpha) %>% as.data.frame() %>% rownames_to_column(var='base') %>% as_data_frame()

    #write tables of significantly DE bases
    write.table(resdf %>% filter(padj < alpha) %>% arrange(padj),
            outtable.spikenorm,
            quote=FALSE,
            sep = "\t",
            row.names=FALSE)

    write.table(resdf.nospike %>% filter(padj < alpha) %>% arrange(padj),
            outtable.libsizenorm,
            quote=FALSE,
            sep = "\t",
            row.names=FALSE)

    #return DESeqDataSet objects
    result = list(dds = dds, dds.nospike = dds.nospike)
    return(result)
}

base = get_de_bases(snakemake@params[["alpha"]], snakemake@params[["samples"]], snakemake@params[["samplegroups"]], snakemake@input[["exp"]], snakemake@input[["si"]], snakemake@output[["spikenorm"]], snakemake@output[["libsizenorm"]])

save(list = ls(all.names = TRUE), file="diff_exp/de_bases/.baseRData", envir = .GlobalEnv)
