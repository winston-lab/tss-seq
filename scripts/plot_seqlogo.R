library(tidyverse)
library(forcats)
library(argparse)
library(ggseqlogo)
#we don't use ggseqlogo plotting because it doesn't allow prior to be taken into account
#we just use their dataframes with font information for geom_polygon

parser = ArgumentParser()
parser$add_argument('-i', '--input', type='character', nargs="+")
parser$add_argument('-t', '--tss_classes', type='character', nargs="+")
parser$add_argument('-b', '--slop', type='integer')
parser$add_argument('-g', '--gc_pct', type='double')
parser$add_argument('-l', '--group_label', type='character', nargs="+")
parser$add_argument('-o', '--out_path', type='character')
args = parser$parse_args()

#from ggseqlogo by Omar Wagih
#https://github.com/omarwagih/ggseqlogo
get_font = function(font){

    GGSEQLOGO_FONT_BASE = getOption('GGSEQLOGO_FONT_BASE')
    if(is.null(GGSEQLOGO_FONT_BASE)){
        GGSEQLOGO_FONT_BASE=system.file("extdata", "", package = "ggseqlogo")
        options(GGSEQLOGO_FONT_BASE=GGSEQLOGO_FONT_BASE)
    }

    #all_fonts = c('sf_bold', 'sf_regular', 'ms_bold', 'ms_regular', 'xkcd_regular')
    font = match.arg(tolower(font), list_fonts(F))
    font_filename = paste0(font, '.font')
    font_obj_name = sprintf('.ggseqlogo_font_%s', font)

    font_obj = getOption(font_obj_name)
    if(is.null(font_obj)){
        # Not loaded into global env yet - load it into options
        font_path = file.path(GGSEQLOGO_FONT_BASE, font_filename)
        font_obj_list = list( tmp=readRDS(font_path) )
        names(font_obj_list) = font_obj_name
        options(font_obj_list)
        font_obj = font_obj_list[[1]]
    }
    # Return font data
    as_tibble(font_obj)
}

main = function(data_paths, tss_classes, slop, gc_pct, group_name, out_path){
    #theoretical maximum (just Shannon entropy)
    max_information = -((1-gc_pct)*log2((1-gc_pct)/2)+gc_pct*log2(gc_pct/2))

    df = tibble()
    for (i in 1:length(data_paths)){
        df = read_tsv(data_paths[i],
                      comment="#",
                      col_names=c('position','A','C','G','T','entropy','low','high','weight')) %>%
            mutate(position = position - as.integer(slop+1)) %>%
            gather(key=base, value=count, c('A','C','G','T')) %>%
            group_by(position) %>%
            mutate(height=entropy*count/sum(count)) %>%
            arrange(height, .by_group=TRUE) %>%
            mutate(base_low = lag(cumsum(height), default=0),
                   base_high = base_low+height) %>%
            left_join(get_font('roboto_bold'), b=c("base"="letter")) %>%
            group_by(position, base) %>%
            mutate(x = scales::rescale(x, to=c((1-(first(weight)))/2, 1-(1-first(weight))/2)),
                   x = x+(position-0.5),
                   y = scales::rescale(y, to=c(first(base_low), first(base_high))),
                   tss_class = tss_classes[i]) %>%
            bind_rows(df, .)
    }
    df = df %>% mutate(tss_class = fct_inorder(tss_class, ordered=TRUE))

    logo = ggplot(data = df, aes(x=x, y=y, group=interaction(position, base), fill=base)) +
        geom_hline(yintercept = max_information, linetype="dashed") +
        geom_polygon() +
        scale_fill_manual(values = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839'),
                          breaks = c('A','C','G','T','U')) +
        scale_y_continuous(limits = c(0, max_information*1.05),
                           breaks = scales::pretty_breaks(n=2),
                           expand=c(0,0),
                           name = "bits") +
        scale_x_continuous(labels = function(x)if_else(x==0, "TSS", as.character(x)),
                           name = "distance from TSS (nt)") +
        facet_grid(tss_class~., switch="y") +
        ggtitle(group_name) +
        theme_classic() +
        theme(legend.position = "none",
              text = element_text(size=12, color="black", face="bold"),
              axis.text.x = element_text(size=12, color="black", face="bold"),
              axis.text.y = element_text(size=10, color="black", face="plain"),
              axis.title.x = element_text(size=10, face="plain"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5),
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=12, face="plain"),
              strip.background = element_blank(),
              strip.text = element_text(size=12),
              strip.text.y = element_text(angle=-180, hjust=1),
              strip.placement = "outside")

    ggsave(out_path, plot=logo, width=14, height=16, units="cm", dpi=326)
}

main(data_paths = args$input,
     tss_classes = args$tss_classes,
     slop = args$slop,
     gc_pct= args$gc_pct,
     group_name = paste(args$group_label, collapse=" "),
     out_path = args$out_path)

