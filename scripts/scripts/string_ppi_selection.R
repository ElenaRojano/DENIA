#! /usr/bin/env Rscript

library(optparse)

################################################################
## OPTPARSE
################################################################

option_list <- list(
        make_option(c("-i", "--input_file"), type="character",
                help="STRING Human PPI translated to ENSEMBL genes (combined scores included)."),
        make_option(c("-c1", "--c1_not_expressed"), type="character",
                help="From expression analysis results, first column with not-expressed genes."),
        make_option(c("-c2", "--c2_not_expressed"), type="character",
                help="From expression analysis results, second column with not-expressed genes."),
        make_option(c("-k", "--combined_score"), type="integer", default=0,
                help="Minimum combined score"),
        make_option(c("-o", "--output_path"), type="character", default="results", 
                help="Output path"))
opt <- parse_args(OptionParser(option_list=option_list))

################################################################
## MAIN
################################################################

# Read human genes interactions from STRING
human_PPI <- read.table(opt$input_file, header = FALSE, sep="\t")

# We kept interactions for which both proteins are expressed in at least one datasets:
not_C1 <- read.table(opt$c1_not_expressed, header = FALSE, sep="\t") # 458
not_C1_vector <- not_C1$V1
filtered_PPI1_c1 <- human_PPI[ ! human_PPI$Ensembl_protein1 %in% not_C1_vector, ]


not_C2 <- read.table(opt$c2_not_expressed, header = FALSE, sep="\t") #6356
not_C2_vector <- not_C2$V1

# we eliminate from the table all the columns that have some element that is not expressed (not_C1_vector)
filtered_PPI1_c2 <- filtered_PPI1_c1[ ! filtered_PPI1_c1$Ensembl_protein2 %in% not_C2_vector, ]

#write.table(filtered_PPI1_c2, "PPI_expressed.csv", sep = "\t", row.names = FALSE)

filtered_PPI <- filtered_PPI1_c2[filtered_PPI1_c2$combined_score >= opt$combined_score, ] 
write.table(filtered_PPI, opt$output_path, sep = "\t", row.names = FALSE)
