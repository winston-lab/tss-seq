library(tidyverse)

melt = function(inmatrix, group, sample, binsize, upstream, downstream, outpath){
    raw = read_table2(inmatrix, skip=3, col_names=FALSE)
    names(raw) = seq(ncol(raw))
    
    df = raw %>%
          rownames_to_column(var="index") %>%
          gather(key = variable, value=value, -index, convert=TRUE) %>%
          transmute(group = group, 
                    sample = sample,
                    index = as.numeric(index),
                    position = (as.numeric(variable)*binsize-upstream)/1000,
                    cpm = as.numeric(value))
    
    write.table(df, file=outpath, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    return(df)
}

melt(inmatrix = snakemake@input[["matrix"]],
     group = snakemake@params[["group"]],
     sample = snakemake@params[["sample"]],
     binsize = snakemake@params[["binsize"]],
     upstream = snakemake@params[["upstream"]],
     downstream = snakemake@params[["dnstream"]],
     outpath = snakemake@output[[1]])