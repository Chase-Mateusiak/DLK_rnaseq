library(tidyverse)
library(readxl)
library(here)

fastq_df = read_tsv(here("data/fastq_lookup.txt"), col_names = 'fastqFileName') %>%
  mutate(Sample = as.numeric(str_extract(fastqFileName, "\\d+")),
         read = ifelse(str_detect(fastqFileName, "R1"), 'one', 'two'))

df = read_excel(here("data/RNAseq_studydesign.xlsx")) %>%
  left_join(fastq_df) %>%
  filter(complete.cases(.)) %>%
  pivot_wider(names_from = read, values_from = fastqFileName) %>%
  dplyr::rename(fastq_1 = one,
                fastq_2 = two,
                sample = Sample,
                pretreatment = Pretreatment,
                treatment = Treatment,
                researcher = `Researcher (RNA isolation)`,
                concentration = Concentration,
                rna_260_280 = `260/280`,
                rna_260_230 = `260/230`) %>%
  mutate(strandedness = 'unstranded') %>%
  dplyr::select(sample, fastq_1, fastq_2, strandedness, pretreatment, treatment, researcher, concentration, rna_260_280, rna_260_230)

# write_csv(df, here("data/nf_sample_sheet.csv"))
