library(tidyverse)
library(here)
library(Mus.musculus)

# set ggplot global options
theme_set(theme_minimal())
theme_update(text = element_text(size=15),
             panel.border = element_rect(colour = "black", fill=NA, size=1))

#'
#' filter tables for padj (less than or equal to) and
#'   log2FoldChange (greater than or equal to), select gene, baseMean,
#'   log2FoldChange, padj
#'
#' @param df a dataframe, expected to be a DESeq results table
#' @param padj_thres 1-alpha threshold. note that filter is less than or equal
#' @param abs_lfc_thres abs(log2fc), note that filter is less than or equal
#'
#' @return a filtered and column subsetted table
#'
heatmap_sig_filter = function(df, padj_thres = .05, abs_lfc_thres = .5){

  df %>%
    filter(padj <= padj_thres,
           abs(log2FoldChange) >= abs_lfc_thres) %>%
    dplyr::select(gene, baseMean, log2FoldChange, padj)

}

#'
#' join a results table with columns SYMBOL in the go_term_table to gene
#'   in the results table. Also filter out any NA log2FoldChange values
#'
#' @param res_df a results table
#' @param go_term_table a table of go terms mapped to gene identifiers
#'  (in this case SYMBOL) created with Mus.musculus package
#'
#' @return go_term_table joined to res_df by SYMBOL = gene
join_res_to_go = function(res_df, go_term_table){

  go_term_table %>%
    left_join(res_df, by = c("SYMBOL" = "gene")) %>%
    filter(!is.na(log2FoldChange))

}

#'
#' given a list of go terms, use Mus.musculus annotation object to get
#'   mapping from go term to ensembl and symbol
#'
#' @param go_terms a list of go terms, eg list(actin: "GO:123456")
#'
#' @return a tibble with columns GO, TERM, ENSEMBL, SYMBOL and maybe a couple others
make_go_table = function(go_terms){

  go_term_table = select(
    Mus.musculus,
    keys = go_terms,
    keytype = "GO",
    columns = c("TERM","ENSEMBL", "SYMBOL")) %>%
    distinct(across(c("ENSEMBL", "GO")), .keep_all = TRUE) %>%
    as_tibble()

}

#'
#' given a list of data and go terms, create a table which has the GO terms
#' and the lfc from the data
#'
#' @param path path to a results table
#' @param go_terms list of go terms, eg
#'   go_terms = list(actin_cytoskeleton = "GO:0015629",
#'   microtubules = "GO:0015630")
#' @param padj_thres padj_threshold to filter results table
#' @param abs_lf_thres abs(log2FoldChange) threshold for filter
#'
#' @return a table of go terms with associated effect values. note that effect
#'  effect values are filtered by padj log2fc
#'
associate_go_with_data = function(path, go_table, padj = .05, abs_lfc_thres = .5){

  data_df_list = map(path, read_csv)

  data_df_list_fltr = map(data_df_list, heatmap_sig_filter)


  out = map(data_df_list_fltr, join_res_to_go, go_table)

  out

}

go_terms = list(
  actin_cytoskeleton = "GO:0015629",
  microtubules = "GO:0015630")

data_paths = list(
  effect_of_dlki_noc = here("data/effect_of_dlki_noc.csv"),
  effect_of_dlki = here("data/effect_of_dlki.csv"),
  effect_of_noc = here("data/effect_of_noc.csv")
)

#'
#' create a heatmap. group rows by up/down reg
#'   MULTIPLY EFFECT BY -1 TO SWTICH SIGNS
#'
#' @param term the go id to plot (assumes table has more than 1 go id)
#' @param df dataframe with appropriate columns
#'
#' @return ggplot heatmap
#'
create_term_heatmap = function(go_id, df){

  fltr_df = df %>%
    filter(go == go_id, condition == "effect_of_dlki_noc")

  term = unique(fltr_df$term)

  symbol_order = fltr_df %>%
    group_by(term) %>%
    arrange(desc(log2FoldChange)) %>%
    pull(symbol) %>%
    unique()

  fltr_df %>%
    mutate(symbol = factor(symbol, levels = symbol_order)) %>%
    mutate(log2FoldChange = -1*log2FoldChange) %>%
    ggplot(aes(condition, symbol, fill = log2FoldChange)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "#075AFF",
                         mid = "white",
                         high = "#FF0000") +
    ggtitle(term)

  # +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

go_terms = list(
  actin_cytoskeleton = "GO:0015629",
  microtubules_cytoskeleton = "GO:0015630")

go_table = make_go_table(go_terms)

heatmap_df_list = map(data_paths, associate_go_with_data, go_table)

heatmap_df_list_mut = map(names(heatmap_df_list),
                          ~mutate(heatmap_df_list[[.]][[1]], condition = .))

heatmap_df = do.call('rbind', heatmap_df_list_mut)

colnames(heatmap_df)[c(1,2:6)] = tolower(colnames(heatmap_df)[c(1,2:6)])

heatmap_by_term = map(go_terms, create_term_heatmap, heatmap_df)








