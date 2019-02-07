library(tidyverse)
library(viridis)

main = function(input_path,
                expression_plot_out,
                shift_plot_out,
                output_path){
    input = read_tsv(input_path)

    df_temp = input %>%
        select(chrom, start, end, name, strand, condition_expr, control_expr) %>%
        group_by(chrom, name, strand) %>%
        mutate(max_condition = max(condition_expr),
               max_control = max(control_expr))

    df_condition = df_temp %>%
        filter(condition_expr == max_condition) %>%
        select(-c(end, max_condition, max_control, control_expr))

    df_control = df_temp %>%
        filter(control_expr==max_control) %>%
        select(-c(end, max_condition, max_control, condition_expr))

    df = df_control %>%
        full_join(df_condition,
                  by=c("chrom", "name", "strand"),
                  suffix = c("_control", "_condition")) %>%
        select(chrom, name, strand, start_control, control_expr, start_condition, condition_expr)

    expression_cutoff = quantile(input[["mean_expr"]], 0.9, na.rm=TRUE)

    expression_plot = ggplot(data = input) +
        geom_vline(xintercept = expression_cutoff,
                   color = "grey65") +
        geom_hline(yintercept = expression_cutoff,
                   color = "grey65") +
        geom_abline(slope = 1,
                    intercept = 0,
                    color = "grey65") +
        stat_bin_hex(aes(x = control_expr + 1,
                         y = condition_expr + 1,
                         color = ..count..,
                         fill = ..count..),
                     geom = "point",
                     size = 0.5,
                     binwidth = c(1e-2, 1e-2),
                     alpha = 0.3) +
        scale_x_log10(name = "control expression") +
        scale_y_log10(name = "condition expression") +
        scale_fill_viridis(guide = FALSE) +
        scale_color_viridis(limits = c(NA, 75),
                            oob = scales::squish,
                            guide = FALSE) +
        theme_light()

    ggsave(expression_plot_out,
           plot=expression_plot,
           width=16,
           height=9,
           units="cm")

    df_changed = df %>%
        filter(control_expr > expression_cutoff &
                   condition_expr > expression_cutoff) %>%
        filter(start_control != start_condition) %>%
        mutate(shift = if_else(strand=="+",
                               start_condition-start_control,
                               start_control-start_condition)) %>%
        select(chrom, name, strand, start_control, start_condition,
               shift, control_expr, condition_expr) %>%
        write_tsv(output_path)

    shift_plot = ggplot(data = df_changed,
           aes(x = shift)) +
        geom_histogram(binwidth = 1,
                       fill="#114477") +
        scale_x_continuous(name = "shift in major genic TSS (nt)",
                           labels = function(x) (if_else(x < 0,
                                                        as.character(x),
                                                        paste0("+", x)))) +
        scale_y_continuous(expand = c(0,0)) +
        theme_light() +
        theme(axis.text = element_text(color="black",
                                       size=11))

    ggsave(shift_plot_out,
           plot=shift_plot,
           width=16,
           height=9,
           units="cm")
}

main(input_path = snakemake@input[["results"]],
     expression_plot_out = snakemake@output[["expression_plot"]],
     shift_plot_out = snakemake@output[["shift_plot"]],
     output_path = snakemake@output[["tsv"]])

