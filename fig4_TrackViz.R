library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# import data 
v_peaks <- import.bed("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/V2/callpeak_results/MACS_results_summits.bed")
d_peaks <- import.bed("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/D2/callpeak_results/MACS_results_summits.bed")
p_peaks <- import.bed("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/P2/callpeak_results/MACS_results_summits.bed")

v_bw <- import.bw("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/V2/coverage.bw", 
                   as = "GRanges")
d_bw <- import.bw("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/D2/coverage.bw", 
                   as = "GRanges")
p_bw <- import.bw("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/P2/coverage.bw", 
                   as = "GRanges")


myAxisTrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "mm10", chromosome = "chrX")

gtrack <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, chromosome = "chrX")

library(biomaRt)
martList <- listMarts()
mart = useMart("ENSEMBL_MART_ENSEMBL")
dataList <- listDatasets(mart)
mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

bgrTrack <- BiomartGeneRegionTrack(genome="mm10",
                                   start=142925000,
                                   end=147957000,
                                   biomart = mart,
                                   chromosome = "chrX",
                                   name="ENSEMBL")

plotTracks(list(bgrTrack), chromosome = "chrX", from = 1e6, to = 1.4e6, transcriptAnnotation = "symbol")

Htr2c_v <- DataTrack(v_bw,chromosome="chrX",
                         from=146925000,
                         to=146957000,
                         name="Vehicle",
                         background.title = "#74A089",
                         col.title = "black",
                         col.axis = "black",
                         type=("hist"))

Htr2c_d <- DataTrack(d_bw,chromosome="chrX",
                     from=146925000,
                     to=146957000,
                     name="DFMO",
                     background.title = "#FDDDA0",
                     col.title = "black",
                     col.axis = "black",
                     type=("hist"))

Htr2c_p <- DataTrack(p_bw,chromosome="chrX",
                     from=146925000,
                     to=146957000,
                     name="Putrescine",
                     background.title = "#F8AFA8",
                     col.title = "black",
                     col.axis = "black",
                     type=("hist"))

plotTracks(c(itrack, Htr2c_v, Htr2c_d, Htr2c_p, gtrack, myAxisTrack),
           chromosome="chrX",
           from=146925000,
           to=146957000,
           transcriptAnnotation = "Name",
           ylim=c(0,250), cex=1, scale=0.25,
           type=("hist"))


availableDisplayPars(Htr2c_v)

