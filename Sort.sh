#!/bin/bash 

R_DIR=/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI


picard=/Users/schvartzmanlab/Desktop/picard/picard.jar

for i in $(cat $R_DIR/Sample_List.txt); do
        echo "Running $i"
        cd $R_DIR/$i
      
		
		samtools sort -@ 4 -o sample_aligned_sorted.bam sample_aligned.bam |
		java -jar $picard MarkDuplicates \
                        INPUT=sample_aligned_sorted.bam \
                        METRICS_FILE=marked_dup_metrics.txt \
                        OUTPUT=aligned_no_duplicates.bam \
                        REMOVE_DUPLICATES=true
		
		
done
