#Ranges assignment - Beatriz Akemi Mizoguchi  
###Importing the BED file of dbSNP variants for mouse chromosome 1
`setwd("~/Desktop/BioDataSkillsClass/bds-files/chapter-09-working-with-range-data")`  
`library(BiocInstaller)`  
`biocLite("GenomicFeatures")`  
`biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")`  
`library(TxDb.Mmusculus.UCSC.mm10.ensGene)`  
`txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene`  
`mm_gtf <- import('Mus_musculus.GRCm38.75_chr1.gtf.gz')`  
`library(rtracklayer)`  
`dbsnp137 <- import("mm10_snp137_chr1_trunc.bed")`  

###Using the reduce() function to extract and collapse all overlapping exons
`collapsed_exons <- reduce(exons(txdb), ignore.strand=TRUE)`

###Subsetting the chromosome 1 collapsed exons
`chr1_collapsed_exons <- collapsed_exons[seqnames(collapsed_exons) == "chr1"]`

###Summary of the dataframe
`summary(width(dbsnp137))`
#### There is one variant that is 732 bases long, which is a bit large compared to the others variants, that are 1.138 mean bases long
`dbsnp137[width(dbsnp137) == 0] #take a look at the variants that are zero bases long`
#####Zero bases long variants do not overlap any feature, so they were not included on the file

###resize() function to include the zero base long features
`dbsnp137_resized <- dbsnp137` #creating a new dataframe  
`zw_i <- width(dbsnp137_resized) == 0` #selecting the variants with zero base long to a variable  
`dbsnp137_resized[zw_i] <- resize(dbsnp137_resized[zw_i], width=1)` #using the resize() function to include the zero base long variants to the set of ranges

###Using findOverlaps() to creat an object called Hits, asking the function to ignore strand
`hits <- findOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)`

###Check Hits and parameters
`hits `  
`length(unique(queryHits(hits)))` #number of dbSNP variants on chromosome 1 that overlap exons  
`length(unique(queryHits(hits)))/length(dbsnp137_resized) `#percentage of dbSNP variants of chromosome 1

###Count the number of variants that overlap each exonic region; the arguments are reversed because the counts need to be based on exons, and the function is asked to ignore strand
`var_counts <- countOverlaps(chr1_collapsed_exons, dbsnp137_resized, ignore.strand=TRUE)`

###Append the variant counts to the chromosome 1 exonic regions
`chr1_collapsed_exons$num_vars <- var_counts`

`chr1_collapsed_exons`

###Export file with column with number of variants per exon appended
`write.table(chr1_collapsed_exons, "chr1_collapsed_exons", sep="\t")`