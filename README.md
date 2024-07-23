# intragenic_drop-off
Scripts for the drop-off probability estimation analysis in the article : FACT maintains chromatin architecture and thereby stimulates RNA polymerase II pausing during transcription in vivo

DOI : https://doi.org/10.1016/j.molcel.2024.05.003

The scripts performs the following operations with the bam files as input.

1. generate transcribed bases RLE tracks (filled up fragments from aligned read pairs)
2. generate coverages for the beginnings and gene body regions of the gene (TSS to TSS+3kb and TSS+3kb to TES respectively)
3. estimate the dropoff rate from the above coverages
4. generate figure 3I in the paper (using ggplot2)

The following packages need to be installed for the script to run
 1. rtracklayer 1.58.0
 2. GenomicAlignments 1.34.0
 3. Rsamtools 2.14.0
 4. ggplot2 3.3.6
 5. dplyr 1.0.9
 6. tidyr 1.2.0
 7. ggpubr 0.4.0
 8. magrittr 2.0.3
 9. doParallel 1.0.17
 

### Extra notes :
1. For this analysis, fastq files were processed with cutadapt as :
    cutadapt -m 40 -q 20 -o $outfolder/sample1_r1.fastq -p $outfolder/sample1_r2.fastq sample1_r1.fastq sample1_r2.fastq
2. Trimmed fastq files were mapped with STAR with following parameters :
    --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $output_folder/$outfile --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.02 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 500000
