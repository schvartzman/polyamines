library(ChIPpeakAnno)
library(ChIPseeker)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rGREAT)
library(GenomicRanges)
library(gprofiler2)
library(msigdbr)
library(fgsea)



# annotate peaks, and plot distribution of peak locations
dfmo_sites_anno <- annotatePeak(dfmo_sites, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
plotAnnoPie(dfmo_sites_anno, cex = 0.9) 

count(v_GR_Anno[v_GR_Anno$annotation == "Promoter (<=1kb)", ])

v_sites_anno <- annotatePeak(v_sites, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
plotAnnoPie(v_sites_anno, cex = 0.9) 

p_sites_anno <- annotatePeak(p_sites, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
plotAnnoPie(p_sites_anno, cex = 0.9) 

# annotating peaks up in dfmo (FC > 2)
d_vs_v_res_deseq_df <- (as.data.frame(d_vs_v_res_deseq))
d_vs_v_res_deseq_df$Fold <- d_vs_v_res_deseq_df$Fold*-1
d_vs_v_res_deseq_df_up <- filter(d_vs_v_res_deseq_df, d_vs_v_res_deseq_df$Fold > 2)
d_vs_v_res_deseq_granges_up <- makeGRangesFromDataFrame(d_vs_v_res_deseq_df_up)

d_vs_v_res_deseq_df_down <- filter(d_vs_v_res_deseq_df, d_vs_v_res_deseq_df$Fold < -2)
d_vs_v_res_deseq_granges_down <- makeGRangesFromDataFrame(d_vs_v_res_deseq_df_down)

v_vs_d_up_anno <- annotatePeak(d_vs_v_res_deseq_granges_up, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
plotAnnoPie(v_vs_d_up_anno, cex = 0.9) 
v_vs_d_up_anno_Df <- as.data.frame(v_vs_d_up_anno)

v_vs_d_down_anno <- annotatePeak(d_vs_v_res_deseq_granges_down, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
plotAnnoPie(v_vs_d_down_anno, cex = 0.9) 

# annotate peaks in transcriptional start sites (TSS; most likely to reflect open promoters), 
# or distal intergenic regions (most likely to reflect open enhancers)

#this requires GRanges object
v_GR_Anno <- as.GRanges(v_sites_anno)
dfmo_GR_Anno <- as.GRanges(dfmo_sites_anno)
put_GR_Anno <- as.GRanges(p_sites_anno)
v_vs_d_GR_Anno <- as.GRanges(v_vs_d_sites_anno)
p_v_v_GR_anno <- as.GRanges(v_vs_p_sites_anno)
c_v_v_GR_anno <- as.GRanges(v_vs_c_sites_anno)


#for TSS peaks 
v_TSS <- v_GR_Anno[abs(v_GR_Anno$distanceToTSS) < 500]
dfmo_TSS <- dfmo_GR_Anno[abs(dfmo_GR_Anno$distanceToTSS) < 500]
p_TSS <- put_GR_Anno[abs(put_GR_Anno$distanceToTSS) < 500]
v_vs_d_TSS <- v_vs_d_GR_Anno[abs(v_vs_d_GR_Anno$distanceToTSS) < 500]
v_vs_p_TSS <- p_v_v_GR_anno[abs(p_v_v_GR_anno$distanceToTSS) < 500]
v_vs_c_TSS <- c_v_v_GR_anno[abs(c_v_v_GR_anno$distanceToTSS) < 500 | ]


# for promoter (<=1kb) peaks 
v_vs_d_promoter <- v_vs_d_GR_Anno[v_vs_d_GR_Anno$annotation == "Promoter (<=1kb)", ]
v_vs_d_promoter_df <- (as.data.frame(v_vs_d_promoter))
v_vs_d_promoter_df$Fold <- v_vs_d_promoter_df$Fold*-1
v_vs_d_promoter_df <- arrange(v_vs_d_promoter_df, Fold, decreasing=TRUE)
sum(v_vs_d_promoter_df$Fold > 2)
sum(v_vs_d_promoter_df$Fold < -2)

#for distal intergenic peaks
v_vs_d_Intergenic <- v_vs_d_GR_Anno[v_vs_d_GR_Anno$annotation == "Distal Intergenic",]
v_vs_d_Intergenic_df <- (as.data.frame(v_vs_d_Intergenic))
v_vs_d_Intergenic_df$Fold <- v_vs_d_Intergenic_df$Fold*-1
v_vs_d_Intergenic_df <- arrange(v_vs_d_Intergenic_df, Fold, decreasing=TRUE)

dfmo_intergenic <- dfmo_GR_Anno[dfmo_GR_Anno$annotation == "Distal Intergenic",]
dfmo_intergenic_df <- (as.data.frame(dfmo_intergenic))
v_vs_d_Intergenic_df$Fold <- v_vs_d_Intergenic_df$Fold*-1
v_vs_d_Intergenic_df <- arrange(v_vs_d_Intergenic_df, Fold, decreasing=TRUE)

sum(v_vs_d_Intergenic_df$Fold > 2)
sum(v_vs_d_Intergenic_df$Fold < -2)

#pull out fold changes and gene names, if that's of interest (example here is shared DI peaks between melanoma and regen)
v_vs_d_TSS_df <- (as.data.frame(v_vs_d_TSS))
v_vs_d_TSS_df$Fold <- v_vs_d_TSS_df$Fold*-1
v_vs_d_TSS_df <- arrange(v_vs_d_TSS_df, Fold)

v_vs_p_TSS_df <- (as.data.frame(v_vs_p_TSS))
v_vs_p_TSS_df$Fold <- v_vs_p_TSS_df$Fold*-1
v_vs_p_TSS_df <- arrange(v_vs_p_TSS_df, Fold, decreasing=TRUE)

v_vs_c_TSS_df <- (as.data.frame(v_vs_c_TSS))
v_vs_c_TSS_df$Fold <- v_vs_c_TSS_df$Fold*-1
v_vs_c_TSS_df <- arrange(v_vs_c_TSS_df, Fold, decreasing=TRUE)

v_vs_d_all_df <- (as.data.frame(v_vs_d_GR_Anno))
v_vs_d_all_df$Fold <- v_vs_d_all_df$Fold*-1
v_vs_d_all_df <- arrange(v_vs_d_all_df, Fold, decreasing=TRUE)


sum(v_vs_d_TSS_df$Fold > 2)
sum(v_vs_d_TSS_df$Fold < -2)

# write tss and distal intergenic peaks as csv 
write.csv(v_vs_d_Intergenic_df, "v_vs_d_distal_intergenic.csv")
write.csv(v_vs_d_TSS_df, "v_vs_d_tss_ecoli_norm_new.csv")

write.csv(v_vs_d_promoter_df, "v_vs_d_promoter.csv")


# prep for gsea - extract fold change and gene symbol, sort
v_vs_d_TSS_gsea <- v_vs_d_TSS_df$Fold
names(v_vs_d_TSS_gsea) <- v_vs_d_TSS_df$SYMBOL
v_vs_d_TSS_gsea <- na.omit(v_vs_d_TSS_gsea)
v_vs_d_TSS_gsea = sort(v_vs_d_TSS_gsea, decreasing=TRUE)

v_vs_p_TSS_gsea <- v_vs_p_TSS_df$Fold
names(v_vs_p_TSS_gsea) <- v_vs_p_TSS_df$SYMBOL
v_vs_p_TSS_gsea <- na.omit(v_vs_p_TSS_gsea)
v_vs_p_TSS_gsea = sort(v_vs_p_TSS_gsea, decreasing=TRUE)

v_vs_c_TSS_gsea <- v_vs_c_TSS_df$Fold
names(v_vs_c_TSS_gsea) <- v_vs_c_TSS_df$SYMBOL
v_vs_c_TSS_gsea <- na.omit(v_vs_c_TSS_gsea)
v_vs_c_TSS_gsea = sort(v_vs_c_TSS_gsea, decreasing=TRUE)

v_vs_d_intergenic_gsea <- v_vs_d_Intergenic_df$Fold
names(v_vs_d_intergenic_gsea) <- v_vs_d_Intergenic_df$SYMBOL
v_vs_d_intergenic_gsea <- na.omit(v_vs_d_intergenic_gsea)
v_vs_d_intergenic_gsea = sort(v_vs_d_intergenic_gsea, decreasing=TRUE)

v_vs_d_all_gsea <- v_vs_d_all_df$Fold
names(v_vs_d_all_gsea) <- v_vs_d_all_df$SYMBOL
v_vs_d_all_gsea <- na.omit(v_vs_d_all_gsea)
v_vs_d_all_gsea = sort(v_vs_d_all_gsea, decreasing=TRUE)

v_vs_d_promoter_gsea <- v_vs_d_promoter_df$Fold
names(v_vs_d_promoter_gsea) <- v_vs_d_promoter_df$SYMBOL
v_vs_d_promoter_gsea <- na.omit(v_vs_d_promoter_gsea)
v_vs_d_promoter_gsea = sort(v_vs_d_promoter_gsea, decreasing=TRUE)

# GSEA 
v_vs_d_gsea_result <- gseGO(v_vs_d_TSS_gsea, ont="BP", keyType = "SYMBOL",
                               minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
                               OrgDb = org.Mm.eg.db, pAdjustMethod = "none")

v_vs_p_gsea_result <- gseGO(v_vs_p_TSS_gsea, ont="BP", keyType = "SYMBOL",
                            minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
                            OrgDb = org.Mm.eg.db, pAdjustMethod = "none")

v_vs_c_gsea_result <- gseGO(v_vs_c_TSS_gsea, ont="BP", keyType = "SYMBOL",
                            minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
                            OrgDb = org.Mm.eg.db, pAdjustMethod = "none")

v_vs_d_intergenic_gsea_result <- gseGO(v_vs_d_intergenic_gsea, ont="BP", keyType = "SYMBOL",
                            minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
                            OrgDb = org.Mm.eg.db, pAdjustMethod = "none")

v_vs_d_all_gsea_result <- gseGO(v_vs_d_all_gsea, ont="BP", keyType = "SYMBOL",
                                       minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
                                       OrgDb = org.Mm.eg.db, pAdjustMethod = "none")

v_vs_d_promoter_gsea_result <- gseGO(v_vs_d_promoter_gsea, ont="CC", keyType = "SYMBOL",
                                minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
                                OrgDb = org.Mm.eg.db, pAdjustMethod = "none", eps = 0)

test_df <- v_vs_d_TSS_df %>% 
  group_by(SYMBOL) %>% 
  summarise(Fold = max(Fold))

#repeat with peaksets of interest!
result_df <- v_vs_d_gsea_result@result

gsea_barplot <- data.frame(v_vs_d_gsea_result$Description, v_vs_d_gsea_result$NES, v_vs_d_gsea_result$pvalue, v_vs_d_gsea_result$p.adjust)
gsea_barplot_plot <- gsea_barplot[order(v_vs_d_gsea_result$pvalue, decreasing=FALSE),]
gsea_barplot_plot <- gsea_barplot_plot[1:45,]


gsea_barplot_plot %>% mutate(name = v_vs_d_gsea_result.Description, desc(v_vs_d_gsea_result.NES)) %>%
  ggplot(aes(v_vs_d_gsea_result.NES, reorder(v_vs_d_gsea_result.Description, v_vs_d_gsea_result.NES), fill = v_vs_d_gsea_result.pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(name = "p value", colours = rev(pals::viridis(1000))) +
  theme_minimal() +
  labs(x = "NES",
       y = "GO Term") +
  #theme_minimal() +
  theme(axis.title.x = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"), 
        plot.title = element_text(size = 20, face = "bold"))+
  ggtitle("GSEA: DFMO vs Vehicle (TSS)")


v_vs_d_gsea_result_DF <- as.data.frame(v_vs_d_gsea_result)
