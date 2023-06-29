library(DESeq2)
library(here)
library(patchwork)


dds = readRDS(here("data/dlk_complete_dds.rds"))

vsd = varianceStabilizingTransformation(dds)

pcaData <- plotPCA(vsd,
                   intgroup=c("pretreatment",
                              "treatment",
                              "researcher"),
                   returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

by_researcher = ggplot(pcaData, aes(PC1, PC2, color=pretreatment, shape=treatment)) +
  geom_point(size=7) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_color_manual(values = c("#6B1AA1", "#A13D03"))+
  theme_bw(base_size = 20) +
  theme(strip.background =element_rect(fill="black"),
        strip.text = element_text(colour = 'white'))+
  facet_wrap(~researcher)

all_pca = ggplot(pcaData, aes(PC1, PC2, color=pretreatment, shape=treatment)) +
  geom_point(size=7) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_color_manual(values = c("#6B1AA1", "#A13D03"))+
  theme_bw(base_size = 20) +
  theme(strip.background =element_rect(fill="black"),
        strip.text = element_text(colour = 'white'))

ptchwrk = (all_pca / by_researcher) + plot_layout(widths = c(15,15), heights = c(7,7))

# ggsave(here("plots/rep_agreement_pca_20220815.png"),
#        device = 'png',
#        width = 25, height = 10)
