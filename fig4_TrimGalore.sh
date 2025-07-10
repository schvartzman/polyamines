#!/bin/bash 

R_DIR=/Users/schvartzmanlab/Desktop/230721_10T_VDPC_ATAC/C1/ngs/release/230711_JUANMA_MAYA_1_MOUSE_LIBRARY_SELF_1000M_PE75_AVITI
R1="_R1.fastq.gz"
R2="_R2.fastq.gz" 


for i in $(cat $R_DIR/Sample_List.txt); do
        echo "Running $i"
        cd $R_DIR/$i
        pwd
		
		trim_galore --paired --cores 4 --quality 15 --fastqc $i$R1 $i$R2
		
		
done
