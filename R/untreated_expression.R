library(DESeq2)
library(tidyverse)

dds = readRDS(here("data/dlk_complete_dds.rds"))
dds = estimateSizeFactors(dds)

mid_expression_filter <- rowSums(edgeR::cpm(counts(dds, normalized = TRUE)) > 3 ) >= 4

dds_fltr = dds[mid_expression_filter, ]

dds_fltr = DESeq(dds_fltr)

goi =   c("kif2a", "kif2b", "kif2c", "kif24",
          "kif18a", "kif18b", "kif19", "kif26b",
          "kif5a", "sarm1", "map3k13", "stmn4",
          "stmn2","nmnat2")
# "kif2b"  "kif2c"  "kif18a" "kif18b" are filtered by the 'mid_expr_fltr' above
# as low expression

untreated_filter_set = function(genes_of_interest, dds_obj){

  gene_fltr = tolower(rowData(dds_obj)$name) %in% genes_of_interest

  subset_dds = dds_obj[gene_fltr,dds_obj$pretreatment == "DMSO" &
                            dds_obj$treatment == "DMSO"]

  mat = counts(subset_dds, normalized = TRUE)

  rownames(mat) = rowData(subset_dds)$name

  mat
}

unfltr_mat = untreated_filter_set(goi, dds)

fltr_mat = untreated_filter_set(goi, dds_fltr)
