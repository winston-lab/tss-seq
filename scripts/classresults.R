#!/usr/bin/Rscript

library(tidyverse) 

args = commandArgs(trailingOnly=TRUE)

classresults = read_tsv(file("stdin"), col_names=FALSE) 
fullresults = read_tsv(args[1]) %>% select(-lfcSE,-stat,-logpval)

out = fullresults %>% inner_join(classresults, by=c("name"="X1")) %>%
        separate(name, into=c('chrom','strand','start','end'), sep="-") %>% 
        mutate_at(vars(strand), funs(if_else(.=="minus", "-", "+")))

out %>% format_tsv(col_names=FALSE) %>% cat()