library(DESeq2)
library(tidyverse)
library(tximport)
library(here)
library(gprofiler2)
library(org.Mm.eg.db)

WRITE = TRUE

files = Sys.glob("/mnt/scratch/dlk_mus/results/star_salmon/*/quant.sf")
names(files) = basename(dirname(files))
file.exists(files)

sample_df = read_csv(here("data/nf_sample_sheet.csv")) %>%
  dplyr::rename(names = sample) %>%
  mutate(names = as.character(names)) %>%
  left_join(tibble(names = names(files), files = files)) %>%
  dplyr::select(-c(fastq_1,fastq_2,strandedness))

# indexDir <- file.path(dir, "Dm.BDGP6.22.98_salmon-0.14.1")
# fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
#               "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")
# gtfPath <- file.path(dir,"Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
#
# makeLinkedTxome(indexDir=indexDir,
#                 source="LocalEnsembl",
#                 organism="Drosophila melanogaster",
#                 release="98",
#                 genome="BDGP6.22",
#                 fasta=fastaFTP,
#                 gtf=gtfPath,
#                 write=FALSE)
#
# se = tximeta(sample_df)

txi = tximport(files,
               type = 'salmon',
               tx2gene = read_tsv("/mnt/scratch/dlk_mus/results/star_salmon/salmon_tx2gene.tsv"))

txi$abundance = txi$abundance[ ,order(match(colnames(txi$abundance), sample_df$names))]
txi$counts = txi$counts[ ,order(match(colnames(txi$counts), sample_df$names))]
txi$length = txi$length [ ,order(match(colnames(txi$length ), sample_df$names))]

stopifnot(identical(sample_df$names, colnames(txi$abundance)))
stopifnot(identical(sample_df$names, colnames(txi$counts)))
stopifnot(identical(sample_df$names, colnames(txi$length)))

sample_df = sample_df %>%
  mutate(researcher = factor(researcher),
         pretreatment = factor(pretreatment, levels = c("DMSO", "DLKi")),
         treatment = factor(treatment, levels = c("DMSO", "Noc")))

dds = DESeqDataSetFromTximport(
  txi = txi,
  colData = sample_df,
  design = ~researcher + pretreatment*treatment
)

id_to_symbol = select(
  org.Mm.eg.db,
  keys = str_remove(rownames(dds), "\\.\\d+$"),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "GENENAME")
) %>%
  as_tibble() %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  dplyr::rename(id = ENSEMBL,
                symbol = GENENAME)

id_convert = gprofiler2::gconvert(str_remove(rownames(dds), "\\.\\d+$"),
                                  organism = 'mmusculus',
                                  mthreshold = 1)

id_df = tibble(
  id = str_remove(rownames(dds), "\\.\\d+$")
) %>%
  left_join(id_convert, by = c('id' = 'input'))


x = rowData(dds) %>%
  as_tibble(rownames = 'id_with_version') %>%
  mutate(id = str_remove(id_with_version, "\\.\\d+$")) %>%
  left_join(id_df)


rowData(dds) = x %>%
  DataFrame(row.names = x$id_with_version)

if(WRITE){
  write_rds(dds, here("data/dlk_complete_dds.rds"))
}
