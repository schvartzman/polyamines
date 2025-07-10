

#!/bin/bash

for i in $(cat $R_DIR/Sample_List.txt); do
        echo "Running $i"
        cd $R_DIR/$i

      
      samtools index -b -@ 4 aligned_no_duplicates.bam bamCoverage -b aligned_no_duplicates.bam -of bigwig -o coverage.bw

done
