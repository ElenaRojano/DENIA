#! /usr/bin/env Rscript

library(optparse)

################################################################
## OPTPARSE
################################################################

option_list <- list(
        make_option(c("-i", "--input_file"), type="character",
                help="STRING Human PPI."),
        make_option(c("-d", "--geneID_file"), type="character",
                help="GeneID dictionary."),
        make_option(c("-o", "--output_path"), type="character", default="results", 
                help="Output path"))

opt <- parse_args(OptionParser(option_list=option_list))

################################################################
## MAIN
################################################################

network <- read.table(opt$input_file, header = TRUE, sep="\t")
eq <- read.table(opt$geneID_file, header = TRUE, sep="\t")

Ensembl_protein1 <- rep(NA, 11511400)
Ensembl_protein2 <- rep(NA, 11511400)

network <- data.frame(network, Ensembl_protein1)
network$Ensembl_protein1 <- eq$GeneID[match(network$protein1, eq$Gene)]
network <- data.frame(network, Ensembl_protein2)
network$Ensembl_protein2 <- eq$GeneID[match(network$protein2, eq$Gene)]

# escribimos una nueva tabla con los resultados
network_Ensembl <- network[,c(4:5,3)]
write.table(network_Ensembl, opt$output_path, sep ="\t", row.names = FALSE)