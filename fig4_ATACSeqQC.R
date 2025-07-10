
library(ATACseqQC)

bamfile <- c("/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/C2/aligned_no_duplicates.bam")

fragSize <- fragSizeDist(bamfile, c("C3"))
estimateLibComplexity(readsDupFreq(bamfile))
