library(biomaRt)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(vsn)
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(apeglm)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(tibble)
library(mebfunctions)
library(fgsea)
library(viridis)
library(ggpubr)
library(ggbreak)
library(ggforce)
library(ggh4x)
library(rstatix)

# colors

Royal2 = c("#9A8822", "#F5CDB4", "#F8AFA8", "#FDDDA0", "#74A089")
Royal1 = c("#899DA4", "#C93312", "#FAEFD1", "#DC863B")

Redon = list(c("#5b859e", "#1e395f", "#75884b", "#1e5a46", "#df8d71", "#af4f2f", "#d48f90", "#732f30", "#ab84a5", "#59385c", "#d8b847", "#b38711"))
colors <- c("darkgrey", "#1e5a46", "#5b859e")


# 10T sgRosa cells were treated with vehicle, 1mM DFMO, or 1mM putrescine for 72h. 
# RNA was extracted for RNAseq.

# load count matrices

counts <- read.table("est_counts_genes_kallisto.txt")


#filter genes that have very low (>10) counts 
counts <- counts %>%
  dplyr::filter(rowSums(.) >= 10)

v_cols <- c(1:3, 13:15)
v_counts <- counts[,v_cols]

no_osteo_cols <- c(1:9, 10:12)
no_osteo_counts <- counts[,no_osteo_cols]

osteo_cols <- c(1:3, 10:12)
osteo_counts <- counts[,osteo_cols]

# create metadata files
samples <- c("JM001", "JM002", "JM003", "JM004", "JM005", "JM006", "JM007", "JM008", "JM009","JM013", "JM014", "JM015", "JM016", "JM017", "JM018", "JM019", "JM020", "JM021")
conditions <- c("Vehicle", "Vehicle", "Vehicle", "DFMO", "DFMO", "DFMO", "Putrescine", "Putrescine", "Putrescine",
                        "Vehicle", "Vehicle", "Vehicle", "DFMO", "DFMO", "DFMO", "Putrescine", "Putrescine", "Putrescine")
genotypes <- c("sgRosa", "sgRosa", "sgRosa", "sgRosa", "sgRosa", "sgRosa","sgRosa", "sgRosa", "sgRosa", 
                       "sgODC1", "sgODC1", "sgODC1", "sgODC1", "sgODC1", "sgODC1", "sgODC1", "sgODC1", "sgODC1")
pa_content <- c("High", "High", "High", "Low", "Low", "Low","High", "High", "High", 
                         "Low", "Low", "Low", "Low", "Low", "Low", "High", "High", "High")


metadata <- data.frame(sample, condition, genotype, pa_content)
rownames(metadata) <- metadata$samples
metadata <- as.matrix(metadata)


# Run DEseq2
dds <- DESeqDataSetFromMatrix(countData=counts, 
                                  colData=metadata, 
                                  design= ~pa_content)

# vst transform, make heatmap 
vsd <- vst(dds, blind = FALSE)
head(assay(dds), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$genotype, vsd$condition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


# pca (fig2)
pby_colors <- c("#1e344e", "#29486c", "#5f7d97", "#a99c5f", "#f4e088", "#f7e9aa", "#643f64", "#8e578e", "#ae87ae")
pby3 <- c("#5f7d97", "#dfbf78", "#8e578e")
pby4 <- c("#5f7d97", "#dfbf78", "#8e578e", "#1e344e")
pby2 <- c("#1e344e", "#c39a4c")
pby2b <- c("#1e344e", "#a3874b")

plotPCA(vsd, intgroup = c("condition", "genotype"))
pcaData <- plotPCA(vsd, intgroup = c( "genotype", "ondition"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$genotype <- factor(genotype, levels = c("sgRosa", "sgODC1"))
pcaData$condition <- factor(condition, levels = c("Vehicle", "DFMO", "Putrescine"))

pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, 
                                shape = genotype)) +
  geom_point(size =3) +
  scale_color_manual(values = pby3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("") +
  #coord_fixed() +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey10", color = "NA"),
        strip.placement = "outside",
        strip.text.x = element_text(color = "white", face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        panel.background = element_blank(),
        panel.grid.major = element_line(size = .25), 
        legend.position = c(0.4, 0.8),
        legend.background = element_rect(color = "grey", fill = "NA"))
pca_plot

save(pca_plot, file = "10T_VDP_RNAseq_pca_final.rdata")
ggsave(plot = pca_plot, "10T_VDP_RNAseq_pca_final.pdf", height = 3, width = 4)

#specify the order of condition factor
dds$condition <- factor(dds$condition, levels = c("Vehicle","DFMO", "Putrescine"))
dds$pa_content <- factor(dds$pa_content, levels = c("Low", "High"))

#run DESeq
dds <- DESeq(dds)
dds_results <- results(dds)
resultsNames(dds)

# low vs. high polyamines 
pa_low_vs_high <- results(dds, contrast=c("pa_content", "Low", "High"))

#add shrunken log fold change estimates 
results <- lfcShrink(
  dds, 
  coef = 2, 
  res = results
)

head(results)
summary(results)

# filter, sort, coerce into df
pa_low_vs_high_df <- pa_low_vs_high %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  dplyr::arrange(log2FoldChange)


# test by plotting 
plotCounts(dds, gene = "Odc1", intgroup = "pa_content")


# prep for gsea
pa_low_vs_high_gsea <- pa_low_vs_high_df$log2FoldChange
names(pa_low_vs_high_gsea) <- pa_low_vs_high_df$gene
pa_low_vs_high_gsea <- na.omit(pa_low_vs_high_gsea)
pa_low_vs_high_gsea = sort(pa_low_vs_high_gsea, decreasing=TRUE)

# run gsea
pa_low_vs_high_gsea_result <- gseGO(pa_low_vs_high_gsea, ont="BP", keyType = "SYMBOL",
                                             minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
                                             OrgDb = org.Mm.eg.db, pAdjustMethod = "none")


# plot 
gsea_barplot <- data.frame(pa_low_vs_high_gsea_result$Description, pa_low_vs_high_gsea_result$NES, pa_low_vs_high_gsea_result$pvalue, pa_low_vs_high_gsea_result$p.adjust)
gsea_barplot_plot <- gsea_barplot[order(pa_low_vs_high_gsea_result$NES, decreasing=TRUE),]
gsea_barplot_plot <- gsea_barplot_plot[1:25,]

# generate deg lists 

pa_low_vs_high_df_sort <- pa_low_vs_high_df %>% 
  dplyr::select(gene, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene) %>% 
  summarize(stat=mean(log2FoldChange)) %>%
  arrange(desc(stat))


pathways.hallmark <- gmtPathways("mh.all.v2023.1.Mm.symbols.gmt.txt")
go.bp <- gmtPathways("m5.go.bp.v2023.1.Mm.symbols.gmt.txt")
go.cc <- gmtPathways("m5.go.cc.v2023.1.Mm.symbols.gmt.txt")
go.mf <- gmtPathways("m5.go.mf.v2023.1.Mm.symbols.gmt.txt")
pathways.celltypes <- gmtPathways("m8.all.v2023.1.Mm.symbols.gmt.txt")
pathways.pathways <- gmtPathways("m2.cp.v2023.1.Mm.symbols.gmt.txt")

fgseaRes <- fgsea(pathways=go.bp, stats=ranks)

fgseaResTidy_osteo <- fgseaRes %>%
  as_tibble() %>%
  arrange(NES)

# plots for fig. 1 and sfig1 
### for JMS

counts <- plotCounts(Rdds, gene = "Odc1", intgroup = "R_condition", returnData = TRUE)
# t_test_stats <- counts %>%
#   filter(R_condition %in% c("Vehicle", "Osteoblast")) %>%
#   t_test(count ~ R_condition, ref.group = "Vehicle", p.adjust.method = "none") %>%
#   add_y_position()

#blues <- brewer.pal(n = 5, name = "Blues")[c(1,4)]
c10T_osteo <- counts %>% 
  filter(., R_condition %in% c("Vehicle", "Osteoblast")) %>% ggplot(., aes(R_condition, count, fill = R_condition)) + 
  #scale_fill_manual(values = blues) +
  scale_fill_brewer(palette = "Blues") +
  #geom_violin(trim = F) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1, alpha = 1/3) +  
  # stat_pvalue_manual(t_test_stats, label = "p.adj", y.position = "y.position",
  #                    hide.ns = T, size = 2.5, tip.length = 0.005, inherit.aes = FALSE) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("Vehicle", "Osteoblast")),
    label = "p.format",
    hide.ns = TRUE,
    size = 2.5,
    tip.length = 0.005
  ) +
  scale_y_continuous(labels = scales::scientific, limits = c(1, 6000)) +
  #scale_color_manual(values=cols) +
  xlab("day") +
  ylab("counts") +
  scale_x_discrete(labels=c("Vehicle" = "0", "Osteoblast" = "3"))+ 
  #theme(legend.key.size = unit(0.5, 'cm')) +
  labs(title='10T1/2 osteoblast') +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey10", color = "NA"),
        strip.placement = "outside",
        strip.text.x = element_text(color = "white", face = "bold"),
        plot.title = element_text(hjust = 0.5),
        #plot.subtitle = element_markdown(hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        #axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(NA),
        panel.grid.minor.y = element_line(NA),
        panel.grid.major.x = element_line(NA),
        panel.grid.minor.x = element_line(NA),
        legend.position = 'none')
c10T_osteo

save(c10T_osteo, file = "RNAseq_10T12_osteoblast.rdata")





