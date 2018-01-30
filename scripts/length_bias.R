library(tidyverse)
library(ggthemes)
library(gridExtra)

main = function(intable, type, direction, condition, control, n_bins=100,
                out1, out2, out3){
    df = read_tsv(intable)
    names(df)[11:14] = c('start', 'end', 'feature_name', 'distance')
    
    df = df %>% mutate(start_relative_pos = distance/(end-start)) %>% 
        filter(start_relative_pos %>% between(0,1)) %>% 
        mutate(bin = cut(start_relative_pos, breaks = seq(0,1,1/n_bins), labels=FALSE))
    n_peaks = nrow(df)
    
    relative_hist_df = df %>% group_by(bin) %>% count()
    relative_permutation_result = relative_hist_df %>% pull(n) %>% chisq.test(simulate.p.value=TRUE, B=2000)
    relative_pval = signif(relative_permutation_result[["p.value"]], 3)
    
    relative_distance_histogram = ggplot(data = df, aes(x=(bin-0.5)/n_bins)) +
        annotate(geom="rect", xmin=0, xmax=1, ymin=0, ymax=round(n_peaks/100), fill="grey70") +
        geom_vline(xintercept=c(0,1), color="black") +
        geom_histogram(breaks = seq(0,1,0.01), fill="#114477", alpha=0.95, size=0) +
        scale_x_continuous(name=NULL, breaks=seq(0,1,0.25),
                           labels=c(if(type=="intragenic"){"ATG"} else if(type=="antisense"){"sense TSS"},
                                    rep("",3),
                                    if(type=="intragenic"){"stop codon"} else if(type=="antisense"){"CPS"}),
                           expand=c(0,0)) +
        ylab("observed peaks") +
        ggtitle(paste0("position of ", type, " TSSs\n", direction, " in ", condition, " vs. ", control),
                subtitle = bquote(permutation ~ test ~ on ~ chi^2 ~ test ~ "statistic: p=" ~ .(relative_pval))) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.margin = margin(.25, 1.25, .25, .25, "cm"),
              axis.title.x = element_text(size=10, face="plain"),
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=10, face="plain"))
    
    relative_obs_over_expected_barplot = ggplot(data = relative_hist_df, aes(x=(bin-0.5)/n_bins, y=n-n_peaks/n_bins)) +
        geom_vline(xintercept=c(0,1), color="black") +
        geom_col(width=1/n_bins, fill="#114477") +
        scale_x_continuous(name=NULL, breaks=seq(0,1,0.25),
                           labels=c(if(type=="intragenic"){"ATG"} else if(type=="antisense"){"sense TSS"},
                                    rep("",3),
                                    if(type=="intragenic"){"stop codon"} else if(type=="antisense"){"CPS"}),
                           expand=c(0,0)) +
        ylab(expression(bold(observed - expected ~ peaks))) +
        ggtitle(paste0("position of ", type, " TSSs\n", direction, " in ", condition, " vs. ", control),
                subtitle = bquote(permutation ~ test ~ on ~ chi^2 ~ test ~ "statistic: p=" ~ .(relative_pval))) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.margin = margin(.25, 1.25, .25, .25, "cm"),
              axis.title.x = element_text(size=10, face="plain"),
              plot.title = element_text(size=12), 
              plot.subtitle = element_text(size=10, face="plain"))
    
    relative_plots = arrangeGrob(relative_distance_histogram, relative_obs_over_expected_barplot, ncol=2)
    
    reference_df = df %>% group_by(feature_name) %>% slice(1) %>% ungroup() %>%
        transmute(end = end-start)
    n_features = nrow(reference_df)
    
    reference_indices = seq(1, sum(reference_df), length.out = n_peaks) %>%
        round() %>% as.integer()
    
    reference_df = reference_df %>% count(end) %>% 
        mutate(start = lag(end, default=0),
               surviving = n_features-cumsum(lag(n, default=0)))
    
    all_positions = list()
    
    #there has got to be a better way to do this...
    for (i in 1:nrow(reference_df)){
        row = reference_df %>% slice(i)
        all_positions[[i]] = rep(row$start:(row$end-1), rep(row$surviving, row$end-row$start)) 
    }
    
    all_positions = unlist(all_positions) 
    
    reference_positions = all_positions[reference_indices]
    observed_positions = df %>% pull(distance) %>% round() %>% as.integer()
    
    reference_hist = hist(reference_positions, plot=FALSE,
                          breaks=seq(0, max(reference_positions, observed_positions), length.out=400))
    observed_hist = hist(observed_positions, plot=FALSE,
                          breaks=seq(0, max(reference_positions, observed_positions), length.out=400))
    absolute_df = tibble(mids = reference_hist[["mids"]],
                         expected = reference_hist[["counts"]],
                         observed = observed_hist[["counts"]])
    
    absolute_permutation_result = chisq.test(absolute_df[["observed"]],
                                             p=(absolute_df[["expected"]]+1e-10)/
                                                 (sum(absolute_df[["expected"]]+1e-10)),
                                             simulate.p.value = TRUE)
    absolute_pval = signif(absolute_permutation_result$p.value, 3)
    
    absolute_distance_freqpoly = ggplot(data = absolute_df %>% gather(type, count, -mids),
                                        aes(x=mids, y=count, color=type)) +
        geom_vline(xintercept = 0, color="grey65", size=1) +
        geom_step() +
        scale_color_ptol() +
        scale_x_continuous(labels = function(x){if_else(x==0,
                                                        if(type=='intragenic'){"ATG"}
                                                        else if(type=='antisense'){"sense TSS"},
                                                        as.character(x/1000))},
                           name=paste("distance from",
                                      if(type=='intragenic'){"ATG"}
                                      else if (type=='antisense'){"sense TSS"} , "(kb)"),
                           limits = if(max(observed_positions)<4000){c(NA, NA)} else {c(NA, 4000)},
                           expand=c(0,0)) +
        ylab("TSS-seq peaks") +
        ggtitle(paste0("position of ", type, " TSSs\n", direction, " in ", condition, " vs. ", control),
                subtitle = bquote(permutation ~ test ~ on ~ chi^2 ~ test ~ "statistic: p=" ~ .(absolute_pval))) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              legend.text = element_text(size=12),
              legend.title = element_blank(),
              legend.position = c(.95,.95),
              legend.justification = c(1,1),
              legend.background = element_blank(),
              plot.title = element_text(size=12), 
              plot.subtitle = element_text(size=10, face="plain"))
    
    absolute_obs_over_exp_barplot = ggplot(data = absolute_df, aes(x=mids, y=observed-expected)) +
        geom_vline(xintercept = 0, color="grey65", size=1) +
        geom_col(fill="#114477", width = max(reference_positions, observed_positions)/400) +
        scale_x_continuous(labels = function(x){if_else(x==0,
                                                        if(type=='intragenic'){"ATG"}
                                                        else if(type=='antisense'){"sense TSS"},
                                                        as.character(x/1000))},
                           name=paste("distance from",
                                      if(type=='intragenic'){"ATG"}
                                      else if (type=='antisense'){"sense TSS"} , "(kb)"),
                           limits = if(max(observed_positions)<4000){c(NA, NA)} else {c(NA, 4000)},
                           expand=c(0,0)) +
        ylab(expression(bold(observed - expected ~ peaks))) +
        ggtitle(paste0("position of ", type, " TSSs\n", direction, " in ", condition, " vs. ", control),
                subtitle = bquote(permutation ~ test ~ on ~ chi^2 ~ test ~ "statistic: p=" ~ .(absolute_pval))) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.title = element_text(size=12), 
              plot.subtitle = element_text(size=10, face="plain"))
    
    absolute_plots = arrangeGrob(absolute_distance_freqpoly, absolute_obs_over_exp_barplot, ncol=2)
    
    relative_foldchange_plot = ggplot(data = df, aes(x=start_relative_pos, y=log2FoldChange)) +
        geom_vline(xintercept = c(0,1), size=1, color="grey65") +
        geom_point(alpha=0.8, shape=16, stroke=0) +
        ylab(expression(bold(log[2] ~ "fold-change"))) +
        scale_x_continuous(name=NULL, breaks=seq(0,1,0.25),
                           labels=c(if(type=="intragenic"){"ATG"} else if(type=="antisense"){"sense TSS"},
                                    rep("",3),
                                    if(type=="intragenic"){"stop codon"} else if(type=="antisense"){"CPS"})) +
        ggtitle(paste0("position vs. fold-change: ", type, " TSSs\n", direction, " in ", condition, " vs. ", control)) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.title = element_text(size=12))
    
    relative_significance_plot = ggplot(data = df, aes(x=start_relative_pos, y=logpadj)) +
        geom_vline(xintercept = c(0,1), size=1, color="grey65") +
        geom_point(alpha=0.8, shape=16, stroke=0) +
        ylab(expression(bold(-log[10] ~ "q-value"))) +
        scale_x_continuous(name=NULL, breaks=seq(0,1,0.25),
                           labels=c(if(type=="intragenic"){"ATG"} else if(type=="antisense"){"sense TSS"},
                                    rep("",3),
                                    if(type=="intragenic"){"stop codon"} else if(type=="antisense"){"CPS"})) +
        ggtitle(paste0("position vs. significance: ", type, " TSSs\n", direction, " in ", condition, " vs. ", control)) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.title = element_text(size=12))
    
    absolute_foldchange_plot = ggplot(data = df, aes(x=distance, y=log2FoldChange)) +
        geom_vline(xintercept = 0, size=1, color="grey65") +
        geom_point(alpha=0.6, shape=16, stroke=0) +
        ylab(expression(bold(log[2] ~ "fold-change"))) +
        scale_x_continuous(labels = function(x){if_else(x==0,
                                                        if(type=='intragenic'){"ATG"}
                                                        else if(type=='antisense'){"sense TSS"},
                                                        as.character(x/1000))},
                           name=paste("distance from",
                                      if(type=='intragenic'){"ATG"}
                                      else if (type=='antisense'){"sense TSS"} , "(kb)"),
                           limits = if(max(observed_positions)<3000){c(NA, NA)} else {c(NA, 3000)}) +
        ggtitle(paste0("position vs. fold-change: ", type, " TSSs\n", direction, " in ", condition, " vs. ", control)) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.title = element_text(size=12))
    
    absolute_significance_plot = ggplot(data = df, aes(x=distance, y=logpadj)) +
        geom_vline(xintercept = 0, size=1, color="grey65") +
        geom_point(alpha=0.6, shape=16, stroke=0) +
        ylab(expression(bold(-log[10] ~ "q-value"))) +
        scale_x_continuous(labels = function(x){if_else(x==0,
                                                        if(type=='intragenic'){"ATG"}
                                                        else if(type=='antisense'){"sense TSS"},
                                                        as.character(x/1000))},
                           name=paste("distance from",
                                      if(type=='intragenic'){"ATG"}
                                      else if (type=='antisense'){"sense TSS"} , "(kb)"),
                           limits = if(max(observed_positions)<3000){c(NA, NA)} else {c(NA, 3000)}) +
        ggtitle(paste0("position vs. significance: ", type, " TSSs\n", direction, " in ", condition, " vs. ", control)) +
        theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black", face="bold"),
              plot.title = element_text(size=12))
    
    foldchange_significance_plots = arrangeGrob(relative_foldchange_plot, relative_significance_plot,
                                                absolute_foldchange_plot, absolute_significance_plot,
                                                ncol=2)
    
    ggsave(out1, plot=relative_plots, width=22, height=10, units="cm")
    ggsave(out2, plot=absolute_plots, width=22, height=10, units="cm")
    ggsave(out3, plot=foldchange_significance_plots, width=24, height=16, units="cm")
}

main(intable = snakemake@input[["peaks"]],
     type= snakemake@wildcards[["ttype"]],
     direction = snakemake@params[["direction"]],
     condition = snakemake@wildcards[["condition"]],
     control = snakemake@wildcards[["control"]],
     n_bins = snakemake@params[["n_bins"]],
     out1 = snakemake@output[["relative"]],
     out2 = snakemake@output[["absolute"]],
     out3 = snakemake@output[["fc_signif"]])
