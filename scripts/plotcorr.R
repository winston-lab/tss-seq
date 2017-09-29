library(tidyverse)
library(GGally)
library(viridis)

main = function(intable, subset, pcount, samplelist, condition, control, outpath){
    print(paste("subset:", subset))
    print(paste("samplelist:", samplelist))
    print(paste("condition:", condition))
    print(paste("control:", control))

    grouplist = c(condition, control)
    
    df = read_tsv(intable, col_names = c('index', 'sample', 'group', 'signal'))
    
    maxsignal = max(df$signal) + pcount

    if(subset=="TRUE"){
        df = df %>% filter(group %in% grouplist)
    }
    df = df %>% filter(sample %in% samplelist) %>%
            select(-group) %>% spread(sample, signal) 
    #c = df %>% select(-index) %>% cor(use="pairwise.complete.obs") %>% as.data.frame() %>%
    #        rownames_to_column(var="xsample") %>% as_tibble() %>%
    #        gather(key=ysample, value=correlation, -xsample)
    #
    #cplot = ggplot(data = c, aes(x=xsample, y = ysample, fill=correlation)) +
    #            geom_raster() +
    #            scale_fill_viridis(option="inferno")
    
    df = df %>% select(-index)
    plots = list()
    
    placeholder = data.frame(x=c(0,1), y=c(0,1))
    
    #for each row
    for (i in 1:ncol(df)){
        #for each column
        for (j in 1:ncol(df)){
            idx = ncol(df)*(i-1)+j
            #upper right (correlation)
            if (i < j){
                c = cor(df[,i], df[,j], use = "complete.obs")
                plot = ggplot(data = placeholder, aes(x, y)) +
                        geom_blank() +
                        annotate("text", x=0.5, y=0.5, label=sprintf("%.2f",round(c,2)), size=5)
                plots[[idx]] = plot
            }
            #top left to bot right diag (density)
            else if (i == j){
                subdf = df %>% select(i) %>% gather(sample, value)
                plot = ggplot(data = subdf, aes(x=(value+pcount))) +
                        geom_density(aes(y=..scaled..), fill="black", size=0) +
                        scale_y_continuous(breaks=c(0,.5,1)) +
                        scale_x_log10(limit = c(pcount, maxsignal)) +
                        annotate("text", x=sqrt(maxsignal)/2, y=0.5,
                                 label=unique(subdf$sample), size=2) 
                plots[[idx]] = plot
            }
            #bottom left (scatter)
            else {
                subdf = df %>% select(i,j) %>% gather(xsample, xvalue, -1) %>%
                            gather(ysample, yvalue, -c(2:3)) %>%
                            filter(!(xvalue == 0 & yvalue == 0))
                plot = ggplot(data = subdf, aes(x=xvalue+pcount, y=yvalue+pcount)) +
                            #geom_point(size=.5, shape=1, alpha=0.3) +
                            geom_hex(aes(fill=log10(..count..)), bins=50) +
                            scale_fill_viridis(option="inferno") +
                            scale_x_log10(limit = c(pcount, maxsignal)) +
                            scale_y_log10(limit = c(pcount, maxsignal))
                plots[[idx]] = plot
            }
        }
    }
    
    mat = ggmatrix(plots, nrow=ncol(df), ncol=ncol(df),
                   xAxisLabels = names(df), yAxisLabels = names(df), switch="both") +
                    theme_light() +
                    theme(axis.text = element_text(size=6),
                          strip.background = element_blank(),
                          strip.text = element_text(size=6, color="black"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.placement="outside",
                          strip.switch.pad.grid = unit(0, "points"),
                          strip.switch.pad.wrap = unit(0, "points"))
    
    ggsave(outpath, mat, width=1.5+ncol(df)*2, height=ncol(df)*4/3, units="cm")
}    

main(intable = snakemake@input[[1]],
     subset = snakemake@params[["subset"]],
     pcount = snakemake@params[["pcount"]],
     samplelist = snakemake@params[["samplelist"]], 
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     outpath = snakemake@output[[1]])

warnings()
