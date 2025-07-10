
# DiffBind for differential peak calling and downstream enrichment analysis 


library(DiffBind)
library(tidyverse)
library(ChIPpeakAnno)
library(ChIPseeker)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(ggridges)
library(edgeR)
library(csaw)

#### Normalization to e. coli spike in reads

ecoli_reads <- c( 27484, 6917, 35221, 6232, 23311, 28589, 23363, 7403, 14921)


# raw library size 
counts <- dba.peakset(dbObj_dba.count, bRetrieve=TRUE)
counts <- as.data.frame(counts)
counts <- counts[,6:14]
libsize <- colSums(counts)
## calculate size factors:
normfactor <- calcNormFactors(object = counts, method = "TMM")
sizefactors <- normfactor * ecoli_reads / 1000000

# recipricol for deeptools scaling
sizefactors_recip <- 1/sizefactors

# exploring different normalization methods 
dba.plotMA(dbObj_dba.analyze, contrast=1, 
           bNormalized=FALSE, sub="Non-Normalized")
dba.plotMA(dbObj_dba.analyze, method=DBA_DESEQ2, sub="DESeq2:lib:full")

dbs <- dba.report(dbObj_dba.analyze, bDB=TRUE, bGain=TRUE, bLoss=TRUE)
dbs$config$factor <- "normalize"
dbs$class[DBA_ID,] <- colnames(dbs$class)[1] <-  "LIB_Full"
dbs$class[DBA_FACTOR,] <- DBA_NORM_LIB
dbs


# RLE normalization 
dbObj_dba.analyze <- dba.normalize(dbObj_dba.analyze, normalize=DBA_NORM_NATIVE)
dbObj_dba.analyze <- dba.analyze(dbObj_dba.analyze)
dba.plotMA(dbObj_dba.analyze, method=DBA_DESEQ2, sub="DESeq2:RLE:RiP")

# norm to reads in peaks 
dbObj_dba.analyze <- dba.normalize(dbObj_dba.analyze, normalize=DBA_NORM_LIB,
                                   library=DBA_LIBSIZE_PEAKREADS, background=FALSE)
dbObj_dba.analyze <- dba.analyze(dbObj_dba.analyze)
dba.plotMA(dbObj_dba.analyze, method=DBA_DESEQ2, sub="DESeq2:lib:RiP")
dba.show(dbObj_dba.analyze,attributes=c(DBA_ID,DBA_FRIP))


# background normalize 
dbObj_dba.analyze <- dba.normalize(dbObj_dba.analyze, method=DBA_ALL_METHODS,
                                                      normalize=DBA_NORM_NATIVE,
                                                      background=TRUE)
dbObj_dba.analyze <- dba.analyze(dbObj_dba.analyze, method=DBA_ALL_METHODS)
dba.show(dbObj_dba.analyze,bContrasts=TRUE)
par(mfrow=c(2,1))
dba.plotMA(dbObj_dba.analyze, method=DBA_EDGER, sub="edgeR:TMM:background")
dba.plotMA(dbObj_dba.analyze, method=DBA_DESEQ2, sub="DESeq2:RLE:background")
par(mfrow=c(1,1))

# loess fit norm 

dbObj_dba.analyze$config$AnalysisMethod <- DBA_EDGER
dbObj_dba.analyze <- dba.normalize(dbObj_dba.analyze, offsets=TRUE)
dbObj_dba.analyze <- dba.analyze(dbObj_dba.analyze)
dba.plotMA(dbObj_dba.analyze)


# CODE FOR FIGURE - E COLI NORM 
# load contrasts and pathways to peaksets 
samples_spike <- read.csv("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/diffbind_samplesheet_spikein.csv")

dbObj_spike <- dba(sampleSheet = samples_spike)
olap.rate <- dba.overlap(dbObj_spike, mode=DBA_OLAP_RATE)
plot(olap.rate, xlab="Overlapping samples", ylab="Overlapping peaks", type="b")

#set up masks based on condition, use only consensus peaksets for downstream analysis
#in this case, I'm using peaks that show up in 2/3 of samples for each condition 
dbObj_consensus <- dba.peakset(dbObj_spike, consensus = DBA_TREATMENT, minOverlap = 0.66)
dbObj_consensus <- dba(dbObj_consensus, mask = dbObj_consensus$masks$Consensus, minOverlap = 1)
consensus_peaks <- dba.peakset(dbObj_consensus, bRetrieve = TRUE)

#visualize overlap of peaks broken down by condition
cons.ol <- dba.plotVenn(dbObj_consensus, mask = dbObj_consensus$masks$Consensus)

#compute count matrix (this takes a long time)
dbObj_dba.count <- dba.count(dbObj_spike, peaks = consensus_peaks, bParallel = FALSE)

# E COLI NORMALIZED COUNTS SAVED HERE - LOAD THIS AND THEN FOLLOW CODE BELOW! 
save(dbObj_dba.count, file = "ATAC_VDP_spikeinnorm_dbaobj.rdata")
load("ATAC_VDP_spikeinnorm_dbaobj.rdata")

## (https://support.bioconductor.org/p/9147040/)

# vector of mapped ecoli reads in each sample - V1, V2, V3, D1, D2, D3, P1, P2, P3
ecoli_reads <- c( 27484, 6917, 35221, 6232, 23311, 28589, 23363, 7403, 14921)

# normalize using ecoli library sizes
dbObj_spike.count <- dba.normalize(dbObj_dba.count, library = ecoli_reads, normalize = DBA_NORM_LIB)
dbObj_spike.count <- dba.analyze(dbObj_spike.count)

save(dbObj_spike.count, file = "ATAC_VDP_spikeinnorm_dbaobj_count.rdata")

dba.plotMA(dbObj_spike.count, sub="E. Coli Spike-in Normalization", bNormalized = T, contrast = 1, bSmooth=FALSE, dotSize=0.5)

#set up contrasts (by condition)
dbObj_dba.contrast <- dba.contrast(dbObj_spike.count, categories = DBA_TREATMENT)

#perform differential analysis
dbObj_dba.analyze <- dba.analyze(dbObj_dba.contrast)

#use MA plots to visualize spread of differential peaks by contrast
dba.plotMA(dbObj_dba.analyze, contrast = 1)

#capture differentially accessible peaks
dbObj_spike_DA <- dba.report(dbObj_spike.count)

pby_colors <- c("#1e344e", "#29486c", "#5f7d97", "#a99c5f", "#f4e088", "#f7e9aa", "#643f64", "#8e578e", "#ae87ae")
pby3 <- c("#5f7d97", "#dfbf78", "#8e578e")
pby4 <- c("#5f7d97", "#dfbf78", "#8e578e", "#1e344e")
pby2 <- c("#1e344e", "#c39a4c")
pby2b <- c("#1e344e", "#a3874b")
colors_troy <- c("#5f7d97", "#1e344e")

# PCA PLOT FOR FIGURE
#run PCA for all peaks in all samples
plot <- dba.plotPCA(dbObj_spike.count, components=1:2, vColors = c("#5f7d97", "#dfbf78", "#8e578e"), 
                    dotSize = 1)
save(plot, file = "10T_VDP_ATAC_PCA.rdata")


#produce heatmaps for each contrast (change contrast # to generate three plots)
hmap<-colorRampPalette(c("red", "black", "blue")) (n=13)
readscores<- dba.plotHeatmap(dbObj_dba.analyze, contrast =3, correlations = FALSE, 
                             scale="row", colScheme = hmap)
plot(dbObj_spike.count)
dba.show(dbObj_dba.analyze,bContrasts=TRUE)

dba.plotVolcano(dbObj_spike.count)



dba.plotMA(dbObj_dba.analyze,bNormalized=FALSE, th=0, sub="Non-Normalized")

dba.plotMA(dbObj_dba.analyze ,bNormalized=TRUE, sub="Normalized: Library Size")

# extract results file - vehicle vs dfmo
d_vs_v_res_deseq <- dba.report(dbObj_spike.count, method=DBA_DESEQ2, contrast = 1, th=1)
p_vs_v_res_deseq <- dba.report(dbObj_spike.count, method=DBA_DESEQ2, contrast = 2, th=1)

# CODE FOR ACCESIBILITY DENSITY PLOT FIGURE
d_vs_v_res_deseq_df <- as.data.frame(d_vs_v_res_deseq)
d_vs_v_res_deseq_df$Fold <- d_vs_v_res_deseq_df$Fold * -1
d_vs_v_res_deseq_df$treatment <- c("DFMO")
d_vs_v_accessibility <- select(d_vs_v_res_deseq_df, c("Fold", "treatment"))

p_vs_v_res_deseq <- as.data.frame(p_vs_v_res_deseq)
p_vs_v_res_deseq$Fold <- p_vs_v_res_deseq$Fold * -1
p_vs_v_res_deseq$treatment <- c("Putrescine")
p_vs_v_accessibility <- select(p_vs_v_res_deseq, c("Fold", "treatment"))

atac_accessibility <- rbind(p_vs_v_accessibility, d_vs_v_accessibility)
atac_accessibility$treatment <- factor(atac_accessibility$treatment, levels = c("DFMO", "Putrescine"))


save(atac_accessibility, file = "10T_atacseq_accessibility_VDP_dataforplot.rdata")

# plot density plots of fold change between conditoins 
plot2 <- atac_accessibility %>%
#%>% filter(treatment == "DFMO") 
  ggplot(aes(Fold)) +
  #geom_density(aes(y = after_stat(scaled), fill = name),
  #              color = "black", linewidth = 1/2, scale = 2, alpha = 0.75) +
  geom_density(aes( y = after_stat(scaled), stat = "density", 
                    fill = treatment, color = treatment, scale = 2), alpha = 1/2) +
  scale_fill_manual(values = c("#dfbf78", "#8e578e")) +
  scale_color_manual(values = c("#dfbf78", "#8e578e"), guide = "none") +
  #scale_color_manual(values = rep("black", 3)) +
  #facet_grid(rows = vars(factor(genotype, levels = c("sgRosa", "sgODC1")))) +
  #scale_y_discrete(expand = expansion(add = c(0.25, 1.75))) +
  scale_x_continuous(limits = c(-3, 3)) +
  ylim(0, NA)+
  coord_cartesian(clip = "off") +
  labs(fill = "Condition", title = "10T1/2 ATACseq",
       x = "Fold Change in Accessibility vs. Vehicle",
       y = "Normalized Peak Frequency") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey30", color = "NA"),
        strip.text = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        #plot.subtitle = element_markdown(hjust = 0.5),
        axis.line = element_line(lineend = "square"),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.y = element_line(NA),
        panel.grid.major.x = element_line(size = 0.5),
        panel.grid.minor.x = element_line(size = 0.25), 
        legend.position = c(0.2, 0.8),
        legend.background = element_rect(fill = "NA"))

plot2

save(plot2, file = "10T_atacseq_accessibility_VDP.rdata")
ggsave(plot = plot2, "10T_atacseq_accessibility_VDP.pdf", height = 3, width = 4)

# just DFMO for keynote
ggsave(plot = plot2, "10T_atacseq_accessibility_VD.pdf", height = 3, width = 4)
ggsave(plot = plot2, "10T_atacseq_accessibility_VD_small.pdf", height = 2, width = 2)


#find overlap peaks by consensus, make lists of peaks unique to each condition
#as a sanity check you can reference these results against the venn diagram numbers (they are the same)
overlap.Peaks <- dba.overlap(dbObj_consensus, dbObj_consensus$masks$Consensus)
v_sites <- overlap.Peaks$onlyA
dfmo_sites <- overlap.Peaks$onlyB
not_dfmo_sites  <- dba.overlap(dbObj_consensus, dbObj_consensus$masks$Vehicle & dbObj_consensus$masks$Putrescine,
                               mode=DBA_OLAP_RATE)

#export as .bed file for HOMER motif analysis 
export.bed(dfmo_sites,con='dfmo_only_peaks_ecoli_norm.bed')
export.bed(d_vs_v_res_deseq,con='v_vs_d_peaks.bed')
p_sites <- overlap.Peaks$onlyC



# code to explore widening of peaks, etc - this did not make it into a figure 
# plot width of peaks
#get peaks found in all replicates of each condition
all_v_peaks <-  dba.peakset(dbObj_spike.count, 4:6, minOverlap=2, bRetrieve=TRUE)
all_d_peaks <- dba.peakset(dbObj_spike, 4:6, minOverlap=2, bRetrieve = T)
all_p_peaks <- dba.peakset(dbObj_consensus, 3, minOverlap=2, bRetrieve = T)


#coerce to plot
all_v_peaks_df <- as.data.frame(all_v_peaks)
all_v_peaks_df$condition <- c("Vehicle")
all_v_peaks_df_plot <- all_v_peaks_df[,c(4, 7)]

all_d_peaks_df <- as.data.frame(all_d_peaks)
all_d_peaks_df$condition <- c("DFMO")
all_d_peaks_df_plot <- all_d_peaks_df[,c(4, 7)]

all_p_peaks_df <- as.data.frame(all_p_peaks)
all_p_peaks_df$condition <- c("Putrescine")
all_p_peaks_df_plot <- all_p_peaks_df[,c(4, 7)]

#merge data frams and plot 
peaksize <- rbind(all_v_peaks_df_plot, all_d_peaks_df_plot, all_p_peaks_df_plot)
peaksize$condition <- factor(peaksize$condition, levels = c("Vehicle", "DFMO", "Putrescine"))
peaksize$condition <- factor(peaksize$condition, levels = c("Vehicle","DFMO","Putrescine"))

sum(peaksize$condition == "DFMO")

p <- ggplot(filter(peaksize, condition %in% c("Vehicle", "DFMO")),
            aes(width,  after_stat(count), fill = condition, color = condition)) +
  geom_density(adjust = 1.5, color = "black", linewidth = 1/2, alpha = 0.75) +
  theme_classic() +
  scale_fill_manual(values = c("#5f7d97", "#dfbf78", "#8e578e")) +
  #facet_wrap(~ condition, ncol=1) +
  xlim(0, 2500) +
  ylab("Count") + 
  xlab("Peak Width") 
  #ggtitle("Distribution of peak width")
p

p <- peaksize %>% 
  ggplot(mapping = aes(x = condition, y = width, fill = condition)) +
  geom_boxplot(stat = "boxplot") +
  #geom_jitter(position = position_jitter(0.25), color = "black", alpha = 0.25) + 
  scale_fill_manual(values = c("#5f7d97", "#dfbf78", "#8e578e")) + 
  #stat_pvalue_manual(dunn_stat,  y.position = c(55000), label = "p.adj.sci", hide.ns = "FALSE", size = 3, tip.length = 0.01, inherit.aes = FALSE) +
  #stat_summary(fun.data = stat_box_data, geom = "text", position = position_dodge(width = 0.9),
  #             hjust = 0.5, vjust = .9, color = "black", size = 3) +
  theme_bw() +
  theme(legend.position="none") +
  xlab("") +
  ylab("Peak Width") +
  theme_bw()
p

#unique peaks
dfmo_peaks_df <- as.data.frame(dfmo_sites)
dfmo_peaks_df$peak_size <- dfmo_peaks_df$end - dfmo_peaks_df$start
dfmo_peaks_df$condition <- c("DFMO")

veh_peaks_df <- as.data.frame(v_sites)
veh_peaks_df$peak_size <- veh_peaks_df$end - veh_peaks_df$start
veh_peaks_df$condition <- c("Vehicle")

#merge two dataframes to make one plot 
v_vs_d_peaksize <- rbind(veh_peaks_df, dfmo_peaks_df)
v_vs_d_peaksize$condition <- factor(v_vs_d_peaksize$condition, levels = c("Vehicle", "DFMO"))

p <- ggplot(v_vs_d_peaksize, aes(width,  after_stat(scaled), fill = condition, color = condition)) +
  geom_density(adjust = 0.75, color = "black", linewidth = 1/2, alpha = 0.7) +
  theme_classic() +
  scale_fill_manual(values = c("#5f7d97", "#dfbf78")) +
  xlim(0, 1500) +
  ylab("Normalized Density") + 
  xlab("Peak Width") +
  ggtitle("Peaks unique to each conditon - distribution of widths")
p

# count number of differentially accessible peaks that are up vs down in dfmo cells
dfmo_up <- sum(p_vs_v_res_deseq_df$Fold > 1.5)
dfmo_down <- sum(p_vs_v_res_deseq_df$Fold < -1.5)


dfmo_down <- sum(all_d_peaks_df$widt < -1.5)



# plotting profiles 
dbObj_spike.count$config$RunParallel <- TRUE
profiles <- dba.plotProfile(dbObj_spike.count, samples=dbObj_spike.count$mask$Putrescine, merge=NULL)

dba.plotProfile(profiles)



