library(Rsubread)
#buildindex(basename="dmel_index",reference="R:\\Downloads\\dmel-all-chromosome-r6.67.fasta.gz")
#curr_sample <- "WT_3"
trim_dir <- "D:\\RNAseq_Project\\data\\trimmed_fastq\\"

align(index = "dmel_index",
      annot.ext = "D:\\RNAseq_Project\\reference\\dmel-all-r6.67.gtf.gz",
      readfile1 = paste0(trim_dir, curr_sample, "_R1.trimmed.fastq.gz"),
      readfile2 = paste0(trim_dir, curr_sample, "_R2.trimmed.fastq.gz"),
      output_file = paste0(curr_sample, ".bam"),
      isGTF = TRUE,
      type = "rna",
      nthreads = 7,
      sortReadsByCoordinates = TRUE)

counts <- featureCounts(paste0(curr_sample, ".bam"),
              annot.ext = "D:\\RNAseq_Project\\reference\\dmel-all-r6.67.gtf.gz",
              isGTFAnnotationFile = TRUE,
              genome = "R:\\Downloads\\dmel-all-chromosome-r6.67.fasta.gz",
              isPairedEnd = TRUE,
              nthreads = 7,
              useMetaFeatures = FALSE,
              GTF.featureType = "gene"
              )

write.csv(counts$counts, paste0(curr_sample,".csv"))
