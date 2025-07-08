#!/bin/bash

# experiment variables
R_DIR=/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI
R1="*_R1_val_1.fq.gz"
R2="*_R2_val_2.fq.gz"


# must include the prefix of the .bt2 files after the last folder 
genome=/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI/genome/mm10/mm10

# *******************************************************************************

# run bowtie2 (paired end), pipe to samtools to save as .bam file
for i in $(cat $R_DIR/Sample_List.txt); do
        echo "Running $i"
        cd $R_DIR/$i

        bowtie2 -p 4 -x $genome -1 $R1 -2 $R2 | samtools view -h -@ 4 -bS - > sample_aligned.bam

done

