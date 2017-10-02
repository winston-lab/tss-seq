library(tidyverse)
library(GGally)
library(viridis)
library(forcats)


main = function(intable, pcount, samplelist, outpath){
    df = intable %>% read_tsv() %>% gather(key=sample, value=signal, -name) %>%
            filter(sample %in% samplelist)
    df$sample = fct_inorder(df$sample, ordered=TRUE)
    df = df %>% spread(sample, signal) %>% select(-name)
    
    maxsignal = max(df) + pcount
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
                        annotate("text", x=0.5, y=0.5, label=sprintf("%.2f",round(c,2)), size=7)
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
                                 label=unique(subdf$sample), size=3) 
                plots[[idx]] = plot
            }
            #bottom left (scatter)
            else {
                #the filtering is a quick hack to avoid the (0,0) bin taking up
                #all of the colorspace
                subdf = df %>% select(i,j) %>% gather(xsample, xvalue, -1) %>%
                            gather(ysample, yvalue, -c(2:3)) %>%
                            #filter(!(xvalue == 0 & yvalue == 0))
                            filter(!(xvalue < 5*pcount & yvalue < 5*pcount))
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
                    theme(axis.text = element_text(size=9),
                          strip.background = element_blank(),
                          strip.text = element_text(size=9, color="black"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.placement="outside",
                          strip.switch.pad.grid = unit(0, "points"),
                          strip.switch.pad.wrap = unit(0, "points"))
    w = 3+ncol(df)*3
    h = 9/16*w
    ggsave(outpath, mat, width=w, height=h, units="cm", dpi=300, type="cairo")
    print(warnings())
}    

main(intable = snakemake@input[[1]],
     pcount = snakemake@params[["pcount"]],
     samplelist = snakemake@params[["samplelist"]],
     outpath = snakemake@output[[1]])
