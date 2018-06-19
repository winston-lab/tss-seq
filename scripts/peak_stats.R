library(tidyverse)
library(forcats)
library(gridExtra)

theme_default = theme_light() +
    theme(text = element_text(size=12, color="black"),
          axis.text = element_text(color="black"),
          plot.title = element_text(size=12))

import = function(df, path, category){
    read_tsv(path) %>%
        select(1:10) %>%
        distinct() %>%
        mutate(category=category) %>%
        bind_rows(df, .) %>%
        return()
}

import_dist = function(df, path, category) {
    read_tsv(path) %>%
        select(1:10, 17) %>%
        magrittr::set_colnames(., c(colnames(.)[1:10], "dist")) %>%
        distinct() %>%
        mutate(category=category) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(genic_path, intra_path, anti_path, conv_path, div_path, inter_path,
                condition, table_out, size_plot_out, dist_plots_out){

    df = tibble() %>%
        import(genic_path, category="genic") %>%
        import(intra_path, category="intragenic") %>%
        import(anti_path, category="antisense") %>%
        import(conv_path, category="convergent") %>%
        import(div_path, category="divergent") %>%
        import(inter_path, category="intergenic") %>%
        mutate(category = fct_inorder(category, ordered=TRUE),
               peak_size = peak_end-peak_start)

    df %>% group_by(category) %>%
        summarise(n = n(),
                  mean_size = mean(peak_size),
                  sd_size = sd(peak_size),
                  median_size = median(peak_size),
                  pct25_size = quantile(peak_size, probs=0.25),
                  pct75_size = quantile(peak_size, probs=0.75),
                  mean_logqval = mean(peak_logqval),
                  sd_logqval = sd(peak_logqval),
                  median_logqval = median(peak_logqval),
                  pct25_logqval = quantile(peak_logqval, probs=0.25),
                  pct75_logqval = quantile(peak_logqval, probs=0.75)) %>%
        mutate_if(is.numeric, funs(signif(., 3))) %>%
        write_tsv(table_out)

    size_plot = ggplot(data = df,
           aes(x=fct_rev(category), y=peak_size)) +
        geom_violin(bw=5,
                    fill="#114477") +
        geom_boxplot(outlier.size = 0, size=0.2, width=0.15, notch=TRUE) +
        scale_y_continuous(name = "peak size (nt)",
                           breaks = scales::pretty_breaks(n=6),
                           sec.axis = sec_axis(~., breaks = scales::pretty_breaks(n=6))) +
        coord_flip() +
        ggtitle(paste("TSS-seq peaks called in", condition)) +
        theme_default +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(size=12),
              panel.grid.major.y = element_blank())
    ggsave(size_plot_out, plot=size_plot, width=12, height=8, units="cm")


    intra_dist_df = tibble() %>%
        import_dist(intra_path, category="intragenic")

    anti_dist_df = tibble() %>%
        import_dist(anti_path, category="antisense")

    conv_div_dist_df = tibble() %>%
        import_dist(conv_path, category="convergent") %>%
        import_dist(div_path, category="divergent") %>%
        mutate(dist = if_else(category=="divergent", -dist, dist))

    intra_dist_plot = ggplot(data = intra_dist_df, aes(x=dist)) +
        geom_vline(xintercept = 0, color="grey65") +
        geom_density(bw=20, fill="#114477") +
        scale_x_continuous(name = "distance from ATG (nt)") +
        ggtitle(paste("intragenic TSSs called in", condition)) +
        theme_default

    anti_dist_plot = ggplot(data = anti_dist_df, aes(x=dist)) +
        geom_vline(xintercept = 0, color="grey65") +
        geom_density(bw=20, fill="#114477") +
        scale_x_continuous(name = "distance from sense TSS (bp)") +
        ggtitle(paste("antisense TSSs called in", condition)) +
        theme_default

    conv_div_dist_plot = ggplot(data = conv_div_dist_df,
           aes(x=dist, group=category, fill=category, color=category)) +
        geom_vline(xintercept = 0, color="grey65") +
        geom_density(bw=10, alpha=0.7) +
        scale_x_continuous(name = "distance from sense TSS (bp)") +
        scale_fill_manual(values = c("#4477AA", "#CC6677")) +
        scale_color_manual(values = c("#4477AA", "#CC6677")) +
        ggtitle(paste("convergent and divergent TSSs called in", condition)) +
        theme_default +
        theme(legend.position = c(0.99, 0.99),
              legend.justification = c(1,1),
              legend.title = element_blank())

    dist_plots = arrangeGrob(intra_dist_plot,
                             anti_dist_plot,
                             conv_div_dist_plot,
                             ncol=1)
    ggsave(dist_plots_out, plot=dist_plots, width=14, height=16, units="cm")
}

main(genic_path = snakemake@input[["genic"]],
     intra_path = snakemake@input[["intragenic"]],
     anti_path = snakemake@input[["antisense"]],
     conv_path = snakemake@input[["convergent"]],
     div_path = snakemake@input[["divergent"]],
     inter_path = snakemake@input[["intergenic"]],
     condition = snakemake@wildcards[["group"]],
     table_out = snakemake@output[["table"]],
     size_plot_out = snakemake@output[["sizes"]],
     dist_plots_out = snakemake@output[["distances"]])

