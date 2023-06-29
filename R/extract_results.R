library(DESeq2)
library(tidyverse)
library(tximport)
library(here)
library(gprofiler2)

dds = readRDS(here("data/dlk_complete_dds.rds"))

mid_expression_filter <- rowSums(edgeR::cpm(counts(dds)) > 3 ) >= 4

dds_fltr = dds[mid_expression_filter, ]

dds_fltr = DESeq(dds_fltr)

dds_fltr_vst = vst(dds_fltr,
                   nsub = 5000,
                   blind = FALSE)

dds_vst_fltr = vst(dds[mid_expression_filter, ], blind = FALSE)
vst_df_fltr = assays(dds_vst_fltr)[[1]] %>%
  as_tibble(rownames = 'gene') %>%
  pivot_longer(-gene, names_to = 'names', values_to = "vst_expr") %>%
  left_join(sample_df)

vst_df_fltr %>%
  ggplot() +
  geom_density(aes(vst_expr, colour = names))

replace_rownames = function(res, dds){

  rownames(res) = rowData(dds)$name

  res
}



fltr_vst_df = assays(dds_fltr_vst)[[1]] %>%
  as_tibble(rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "vst") %>%
  mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) %>%
  left_join(dplyr::select(id_convert, target, name, description), by=c("gene_id" = "target")) %>%
  mutate(description = str_remove(description, '\\[Source:MGI Symbol;Acc:MGI:95773\\]')) %>%
  dplyr::rename(symbol = name) %>%
  dplyr::rename(names = sample)

researcher_factor = colData(dds_fltr_vst) %>%
  as_tibble() %>%
  pull(researcher) %>%
  factor()

dds_fltr_vst_batch_removed =
  limma::removeBatchEffect(as.matrix(assays(dds_fltr_vst)[[1]]), researcher_factor)

# par(mfrow=c(1,1))
# boxplot(as.data.frame(as.matrix(assays(dds_fltr_vst)[[1]])),main="Original")
# boxplot(as.data.frame(dds_fltr_vst_batch_removed),main="Batch corrected")

batch_removed_fltr_vst_df = dds_fltr_vst_batch_removed %>%
  as_tibble(rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "vst") %>%
  mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) %>%
  left_join(dplyr::select(id_convert, target, name, description), by=c("gene_id" = "target")) %>%
  mutate(description = str_remove(description, '\\[Source:MGI Symbol;Acc:MGI:95773\\]')) %>%
  dplyr::rename(symbol = name) %>%
  dplyr::rename(names = sample) %>%
  left_join(sample_df)


res_list = list(
  effect_of_dlki_noc = results(
    dds_fltr,
    contrast = list( c("pretreatment_DLKi_vs_DMSO", "pretreatmentDLKi.treatmentNoc")),
    test = "Wald"
  ),
  effect_of_noc_dlki = results(
    dds_fltr,
    contrast = list( c("treatment_Noc_vs_DMSO","pretreatmentDLKi.treatmentNoc"))
  ),
  effect_of_dlki = results(
    dds_fltr,
    contrast = c('pretreatment', 'DLKi', 'DMSO')
  ),
  effect_of_noc = results(
    dds_fltr,
    contrast = c('treatment', 'Noc', 'DMSO')
  )
)

shrunken_res_list = map(res_list,
                        ~lfcShrink(dds = dds_fltr, res = ., type = "ashr"))

shrunken_res_list_symbol = map(shrunken_res_list, replace_rownames, dds_fltr)

# shrunken_res_list$effect_of_noc_dlki %>%
#   as_tibble(rownames = 'gene') %>%
#   arrange(log2FoldChange) %>%
#   head() %>%
#   view()
#
#
# vst_df = assays(dds_vst)[[1]] %>%
#   as_tibble(rownames = 'gene') %>%
#   pivot_longer(-gene, names_to = 'names', values_to = "vst_expr") %>%
#   left_join(sample_df)
#
# plt = vst_df %>%
#   filter(gene == 'ENSMUSG00000022044.15') %>%
#   ggplot(aes(pretreatment, vst_expr, color = treatment)) +
#   geom_point()
#
# plotly::ggplotly(plt)
#
# vst_df %>%
#   ggplot() +
#   geom_density(aes(vst_expr, colour = names))

# map(names(shrunken_res_list_symbol),
#     ~write_csv(as_tibble(shrunken_res_list_symbol[[.]], rownames = 'gene'),
#                                            file.path(here("data"),
#                                                      paste0(.,".csv"))))
