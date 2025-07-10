# library(biomaRt)
# library(DESeq2)
# library(EnhancedVolcano)
# library(pheatmap)
# library(vsn)
# library(data.table)
library(tidyverse)
# library(apeglm)
# library(dplyr)
# library(RColorBrewer)
library(ggplot2)
# library(org.Mm.eg.db)
# library(clusterProfiler)
# library(tibble)
library(mebfunctions)
library(fgsea)
library(viridis)
library(ggpubr)
library(ggbreak)
library(ggforce)
library(ggh4x)
library(rstatix)
library(ggtext)
library(sysfonts)
library(extrafont)
library(showtext)

load("/Volumes/chrometab/MEB/RNAseq/230721_10T_VDPB/low_vs_high_pa_gsea_dataframe_forJMS.rdata")
load("/Volumes/chrometab/MEB/RNAseq/230721_10T_VDPB/high_vs_low_pa_gsea_dataframe_forJMS.rdata")
load("/Volumes/chroMetab/MEB/RNAseq/230721_10T_VDPB/low_vs_high_pa_gsea_waterfall.rdata")

###low vs high

fgseaResTidy_low_v_high_mod <- fgseaResTidy_low_v_high %>% arrange(-NES) %>% mutate(pathway = gsub("GOBP_", "", pathway, perl = TRUE)) %>%
  mutate(pathway = gsub(pattern = "GOCC_", replacement = "", pathway, ignore.case = F)) %>%
  mutate(pathway = gsub(pattern = "_", replacement = " ", pathway)) %>%
  mutate(pathway = gsub(pattern = "plus", replacement = "+", pathway)) %>%
  mutate(pathway = tolower(pathway))

fgseaResTidy_low_v_high_mod$pathway_short <- fgseaResTidy_low_v_high_mod$pathway

bold_rows <- c(1,8,10,12,13,15,17,18)
fgseaResTidy_low_v_high_mod <- fgseaResTidy_low_v_high_mod %>% mutate(pathway = ifelse(row_number() %in% bold_rows, paste0("**<span style = 'font-weight:bold; color:black'>", pathway, "</span>**"),
                                                                                       paste0("<span style = 'font-weight:300; color:grey'>", pathway, "</span>"))) %>% 
  mutate(bold = ifelse(row_number() %in% bold_rows, "bold", "italic"))

gsea_plot_low_v_high <- meb_GSEAbarplot(fgseaResTidy_low_v_high_mod, n_terms = 18, fill = T, metric = "NES") +
  ggtitle("GSEA - Low Polyamine vs. High Polyamine") +
  scale_fill_gradientn(name = "p value", colours = rev(pals::viridis(1000))) +
  #scale_y_discrete(labels = y_labels_low_vs_high) +
  #theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        #plot.subtitle = element_markdown(hjust = 0.5),
        #axis.line = element_line(lineend = "square"),
        #axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        panel.background = element_blank(),
        axis.text.y = element_markdown(),
        # panel.grid.major.y = element_line(NA),
        # panel.grid.minor.y = element_line(NA),
        # panel.grid.major.x = element_line(NA),
        # panel.grid.minor.x = element_line(NA),
        #legend.position = 'none'
  )
gsea_plot_low_v_high

ggsave(plot = gsea_plot_low_v_high, "low_vs_high_pa_gsea_barplot_jms.pdf", height = 5, width = 8)

### alt
gsea_plot_low_v_high_alt <- meb_GSEAbarplot(fgseaResTidy_low_v_high_mod, n_terms = 18, fill = T, metric = "NES") +
  ggtitle("GSEA - Low Polyamine vs. High Polyamine") +
  #scale_fill_gradientn(name = "p value", colours = rev(pals::viridis(1000))) +
  scale_fill_viridis(option = "D", direction = -1, begin = 0.4, end = 1) +
  geom_text(aes(label = pathway_short, x = 0.01, fontface = bold),
            position = position_nudge(x = 0.01),
            hjust = 0, # Adjust to ensure text is inside the bar
            color = "black", # Text color
            size = 3.5, # Text size
            family = "Helvetica") + # Text font
  #scale_y_discrete(labels = y_labels_low_vs_high) +
  #theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        #plot.subtitle = element_markdown(hjust = 0.5),
        #axis.line = element_line(lineend = "square"),
        #axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        # panel.grid.major.y = element_line(NA),
        # panel.grid.minor.y = element_line(NA),
        # panel.grid.major.x = element_line(NA),
        # panel.grid.minor.x = element_line(NA),
        #legend.position = 'none'
  )
gsea_plot_low_v_high_alt

ggsave(plot = gsea_plot_low_v_high_alt, "low_vs_high_pa_gsea_barplot_alt_jms.pdf", height = 5, width = 8)

#### high_v_low

fgseaResTidy_high_v_low_mod <- fgseaResTidy_high_v_low %>% mutate(pathway = gsub("GOBP_", "", pathway, perl = TRUE)) %>%
  mutate(pathway = gsub(pattern = "GOCC_", replacement = "", pathway, ignore.case = F)) %>%
  mutate(pathway = gsub(pattern = "_", replacement = " ", pathway)) %>%
  mutate(pathway = gsub(pattern = "plus", replacement = "+", pathway)) %>%
  mutate(pathway = tolower(pathway))

fgseaResTidy_high_v_low_mod$pathway_short <- fgseaResTidy_high_v_low_mod$pathway

gsea_plot_high_v_low <- meb_GSEAbarplot(fgseaResTidy_high_v_low_mod, n_terms = 18, fill = T, metric = "NES") +
  ggtitle("GSEA - High Polyamine vs. Low Polyamine") +
  scale_fill_gradientn(name = "p value", colours = rev(pals::viridis(1000))) +
  #scale_y_discrete(labels = y_labels_low_vs_high) +
  #theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        #plot.subtitle = element_markdown(hjust = 0.5),
        #axis.line = element_line(lineend = "square"),
        #axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        panel.background = element_blank(),
        #axis.text.y = element_markdown(),
        # panel.grid.major.y = element_line(NA),
        # panel.grid.minor.y = element_line(NA),
        # panel.grid.major.x = element_line(NA),
        # panel.grid.minor.x = element_line(NA),
        #legend.position = 'none'
  )
gsea_plot_high_v_low

ggsave(plot = gsea_plot_high_v_low, "high_vs_low_pa_gsea_barplot_jms.pdf", height = 5, width = 8)

gsea_plot_high_v_low_alt <- meb_GSEAbarplot(fgseaResTidy_high_v_low_mod, n_terms = 18, fill = T, metric = "NES") +
  ggtitle("GSEA - High Polyamine vs. Low Polyamine") +
  scale_fill_gradientn(name = "p value", colours = rev(pals::viridis(1000))) +
  scale_fill_viridis(option = "D", direction = -1, begin = 0.4, end = 1) +
  scale_y_discrete(expand = c(0,0)) +
  geom_text(aes(label = pathway_short, x = 0.01),
            position = position_nudge(x = 0.01),
            hjust = 0, # Adjust to ensure text is inside the bar
            color = "black", # Text color
            size = 3.5, # Text size
            family = "Helvetica") + # Text font
  #scale_y_discrete(labels = y_labels_low_vs_high) +
  #theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        #plot.subtitle = element_markdown(hjust = 0.5),
        #axis.line = element_line(lineend = "square"),
        #axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        # panel.grid.major.y = element_line(NA),
        # panel.grid.minor.y = element_line(NA),
        # panel.grid.major.x = element_line(NA),
        # panel.grid.minor.x = element_line(NA),
        #legend.position = 'none'
  )
gsea_plot_high_v_low_alt

ggsave(plot = gsea_plot_high_v_low_alt, "high_vs_low_pa_gsea_barplot_alt_jms.pdf", height = 5, width = 8)

save(gsea_plot_low_v_high,
     gsea_plot_high_v_low,
     gsea_plot_low_v_high_alt,
     gsea_plot_high_v_low_alt,
     file = "high_vs_low_pa_gsea_barplot_jms.rdata")

### waterfall plot

font_add(family = "Helvetica", regular = "/System/Library/Fonts/Helvetica.ttc")
showtext_auto()
loadfonts(device = "pdf")

fgseaResTidy_low_v_high_wfall <- fgseaResTidy_low_v_high_wfall %>%
  mutate(normalized_padj = scales::rescale(-padj, to = c(0.1, 1)))  # Scaling to a range that makes sense for alpha

plot <- filter(fgseaResTidy_low_v_high_wfall, category != "Other") %>% ggplot(aes(reorder(ID, -NES), NES, color = category)) +
  geom_jitter(aes(size = -padj, alpha = normalized_padj), width = 0.05, height = 0.05) +
  #geom_hline(yintercept=c(1.5, -1.5), linetype="dashed", color = "grey") +
  scale_color_manual(values = c("#23334C", "#FFB500", "grey")) +
  scale_alpha_continuous(range = c(1/20, 1/5)) +  # Set alpha limits
  #scale_size_manual(values = c("Cell Cycle" = 5, "Differentiation"= 5, "Other" = 5)) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  guides(alpha = "none", color = guide_legend(override.aes = list(size = 3), order = 1, position = "inside"), size = guide_legend(override.aes = list(alpha = 1), order = 2, position = "inside")) +
  #labs(x = expression("Pathway Rank: \u2190 Low Polyamines | \u2192 High Polyamines"), y = "Normalized Enrichment Score", title = "GSEA category rank distribution") +
  labs(x = paste0("Pathway Rank: ","<span style='font-family:Arial;'>\u2190</span>",
                  " Low Polyamines | High Polyamines ",
                  "<span style='font-family:Arial;'>\u2192</span>"),
       y = "Normalized Enrichment Score",
       title = "GSEA category rank distribution",
       color = "",
       size = "Adjusted p-value") +
  #facet_grid(~ category) +
  #theme_classic() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x=element_blank(),
        axis.title.y = element_text(),
        axis.title.x = element_markdown(),
        axis.ticks.x=element_blank(),
        #plot.subtitle = element_markdown(hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position.inside = c(0.75, 0.7),
        #legend.justification.top = "center",
        #legend.location = "plot",
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = "NA"))
plot
waterfall_plot <- plot


ggsave(plot = plot, "GSEA_category_waterfall.pdf", device = cairo_pdf(family = "Arial"), height = 3, width = 5)
save(waterfall_plot,
     file = "GSEA_waterfall.rdata")

plot_by_category <- fgseaResTidy_low_v_high_wfall %>% ggplot(aes(reorder(ID, -NES), NES, color = category)) +
  geom_jitter(aes(size = -padj, alpha = normalized_padj), width = 0.05, height = 0.05) +
  #geom_hline(yintercept=c(1.5, -1.5), linetype="dashed", color = "grey") +
  scale_color_manual(values = c("#23334C", "#FFB500", "grey")) +
  scale_alpha_continuous(range = c(1/20, 1/5)) +  # Set alpha limits
  #scale_size_manual(values = c("Cell Cycle" = 5, "Differentiation"= 5, "Other" = 5)) +
  scale_x_discrete(expand = c(0.075, 0.075)) +
  guides(alpha = "none", color = guide_legend(override.aes = list(size = 3), order = 1, position = "right"), size = guide_legend(override.aes = list(alpha = 1), order = 2, position = "right")) +
  #labs(x = expression("Pathway Rank: \u2190 Low Polyamines | \u2192 High Polyamines"), y = "Normalized Enrichment Score", title = "GSEA category rank distribution") +
  labs(x = paste0("Pathway Rank: ","<span style='font-family:Arial;'>\u2190</span>",
                  " Low Polyamines | High Polyamines ",
                  "<span style='font-family:Arial;'>\u2192</span>"),
       y = "Normalized Enrichment Score",
       title = "GSEA category rank distribution",
       color = "Category",
       size = "Adjusted p-value") +
  facet_grid2(~ category, strip = strip) +
  #theme_classic() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x=element_blank(),
        axis.title.y = element_text(),
        axis.title.x = element_markdown(),
        axis.ticks.x=element_blank(),
        #plot.subtitle = element_markdown(hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position.inside = c(0.75, 0.7),
        #legend.justification.top = "center",
        #legend.location = "plot",
        legend.box = "vertical",
        legend.background = element_rect(fill = "white", color = "NA"))

plot_by_category

ggsave(plot = plot_by_category, "GSEA_waterfall_bycategory.pdf", height = 4, width = 10)

























