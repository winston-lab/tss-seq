library(tidyverse)
library(broom)
library(viridis)

main = function(condition, control,
                genic_path, intra_path, anti_path, conv_path, div_path,
                tsv_out, lfc_v_lfc_out, lfc_v_expr_out, expr_v_expr_out){
    df = read_tsv(intra_path) %>%
        select(peak_name, mean_expr, log2_foldchange, log10_padj,
               condition_expr, control_expr, feat_name=orf_name,
               dist=atg_to_peak_dist) %>%
        mutate(tclass = "intragenic") %>%
        bind_rows(read_tsv(anti_path) %>%
                      select(peak_name, mean_expr, log2_foldchange, log10_padj,
                             condition_expr, control_expr,
                             feat_name=transcript_name,
                             dist=sense_tss_to_peak_dist) %>%
                      mutate(tclass = "antisense")) %>%
        bind_rows(read_tsv(conv_path) %>%
                      select(peak_name, mean_expr, log2_foldchange, log10_padj,
                             condition_expr, control_expr,
                             feat_name=transcript_name,
                             dist=sense_tss_to_peak_dist) %>%
                      mutate(tclass = "convergent")) %>%
        bind_rows(read_tsv(div_path) %>%
                      select(peak_name, mean_expr, log2_foldchange, log10_padj,
                             condition_expr, control_expr,
                             feat_name=transcript_name,
                             dist=sense_tss_to_peak_dist) %>%
                      mutate(tclass = "divergent")) %>%
        inner_join(read_tsv(genic_path) %>%
                       select(mean_expr, log2_foldchange, log10_padj,
                              condition_expr, control_expr,
                              feat_name=genic_name),
                   by = "feat_name", suffix = c("_class", "_genic")) %>%
        group_by(tclass, feat_name) %>%
        arrange(desc(log10_padj_class), desc(log10_padj_genic), .by_group=TRUE) %>%
        slice(1) %>%
        write_tsv(tsv_out)

    lfc_v_lfc = ggplot() +
        geom_hline(yintercept = 0, size=0.5, color="grey65") +
        geom_vline(xintercept = 0, size=0.5, color="grey65") +
        stat_binhex(data = df,
                    aes(x=log2_foldchange_class, y=log2_foldchange_genic,
                        color=..count..),
                    geom="point", binwidth=c(0.1, 0.1),
                    size=0.5, alpha=0.6, fill=NA) +
        geom_smooth(data = df,
                    aes(x=log2_foldchange_class, y=log2_foldchange_genic),
                    method="lm", size=0.8, color="#114477", alpha=0.8) +
        geom_label(data = df %>%
                       group_by(tclass) %>%
                       do(tidy(lm(.$log2_foldchange_genic ~
                                      .$log2_foldchange_class))) %>%
                       filter(term != "(Intercept)") %>%
                       mutate(label = paste("slope =", signif(estimate, 2),
                                            "\np =",
                                            scales::scientific(p.value))),
                   aes(label=label),
                   x = min(df[["log2_foldchange_class"]]),
                   y = max(df[["log2_foldchange_genic"]]),
                   hjust=0, vjust=1, alpha=0.6, size=10/72*25.4,
                   label.size=0) +
        facet_wrap(~tclass, ncol=2) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        xlab(expression(TSS ~ class ~ log[2] ~ textstyle(frac(condition, control)))) +
        ylab(expression(atop(genic ~ TSS, ~ log[2] ~ textstyle(frac(condition, control))))) +
        ggtitle("TSS-seq: TSS classes vs. genic") +
        theme_light() +
        theme(text = element_text(size=12, color="black"),
              plot.title = element_text(size=12),
              strip.background = element_blank(),
              strip.text = element_text(color="black"),
              axis.text = element_text(color="black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
    ggsave(lfc_v_lfc_out, plot=lfc_v_lfc, width=16, height=12, units="cm")

    lfc_v_expr = ggplot() +
        geom_hline(yintercept = 0, size=0.5, color="grey65") +
        geom_vline(xintercept = 0, size=0.5, color="grey65") +
        stat_binhex(data = df,
                    aes(x=control_expr_genic+1, y=log2_foldchange_class,
                        color=..count..),
                    geom="point", binwidth=c(0.05, 0.1),
                    size=0.5, alpha=0.6, fill=NA) +
        geom_smooth(data = df,
                    aes(x=control_expr_genic+1, y=log2_foldchange_class),
                    method="lm", size=0.8, color="#114477", alpha=0.8) +
        geom_label(data = df %>%
                       group_by(tclass) %>%
                       do(tidy(lm(.$log2_foldchange_class ~
                                      log10(.$control_expr_genic+1)))) %>%
                       filter(term != "(Intercept)") %>%
                       mutate(label = paste("slope =", signif(estimate, 2),
                                            "\np =",
                                            scales::scientific(p.value))),
                   aes(label=label),
                   x = min(df[["control_expr_genic"]]),
                   y = min(df[["log2_foldchange_class"]]),
                   hjust=0, vjust=0, alpha=0.6, size=10/72*25.4,
                   label.size=0) +
        facet_wrap(~tclass, ncol=2) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        scale_x_log10() +
        xlab(paste(control, "genic TSS expression")) +
        ylab(expression(atop(TSS ~ class, ~ log[2] ~ textstyle(frac(condition, control))))) +
        ggtitle("TSS-seq: TSS classes vs. genic") +
        theme_light() +
        theme(text = element_text(size=12, color="black"),
              plot.title = element_text(size=12),
              strip.background = element_blank(),
              strip.text = element_text(color="black"),
              axis.text = element_text(color="black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
    ggsave(lfc_v_expr_out, plot=lfc_v_expr, width=16, height=12, units="cm")

    expr_v_expr = ggplot() +
        geom_hline(yintercept = 0, size=0.5, color="grey65") +
        geom_vline(xintercept = 0, size=0.5, color="grey65") +
        stat_binhex(data = df,
                    aes(x=control_expr_genic+1, y=condition_expr_class+1,
                        color=..count..),
                    geom="point", binwidth=c(0.01, 0.01),
                    size=0.5, alpha=0.6, fill=NA) +
        geom_smooth(data = df,
                    aes(x=control_expr_genic+1, y=condition_expr_class+1),
                    method="lm", size=0.8, color="#114477", alpha=0.8) +
        geom_label(data = df %>%
                       group_by(tclass) %>%
                       do(tidy(lm(log10(.$condition_expr_class+1) ~
                                      log10(.$control_expr_genic+1)))) %>%
                       filter(term != "(Intercept)") %>%
                       mutate(label = paste("slope =", signif(estimate, 2),
                                            "\np =",
                                            scales::scientific(p.value))),
                   aes(label=label),
                   x = min(df[["control_expr_genic"]]),
                   y = log10(max(df[["condition_expr_class"]])),
                   hjust=0, vjust=1, alpha=0.6, size=10/72*25.4,
                   label.size=0) +
        facet_wrap(~tclass, ncol=2) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        scale_x_log10() +
        scale_y_log10() +
        xlab(paste(control, "genic TSS expression")) +
        ylab(paste(condition, "TSS class expression")) +
        ggtitle("TSS-seq: TSS classes vs. genic") +
        theme_light() +
        theme(text = element_text(size=12, color="black"),
              plot.title = element_text(size=12),
              strip.background = element_blank(),
              strip.text = element_text(color="black"),
              axis.text = element_text(color="black"))
    ggsave(expr_v_expr_out, plot=expr_v_expr, width=16, height=12, units="cm")
}

main(condition=snakemake@wildcards[["condition"]],
     control=snakemake@wildcards[["control"]],
     genic_path=snakemake@input[["genic"]],
     intra_path=snakemake@input[["intragenic"]],
     anti_path=snakemake@input[["antisense"]],
     conv_path=snakemake@input[["convergent"]],
     div_path=snakemake@input[["divergent"]],
     tsv_out=snakemake@output[["tsv"]],
     lfc_v_lfc_out=snakemake@output[["lfc_v_lfc"]],
     lfc_v_expr_out=snakemake@output[["lfc_v_expr"]],
     expr_v_expr_out=snakemake@output[["expr_v_expr"]])
