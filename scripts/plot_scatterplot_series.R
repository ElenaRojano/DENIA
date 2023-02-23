#! /usr/bin/env Rscript
# x,y graph

library(ggplot2)
library(optparse)

################################################################
# OPTPARSE
################################################################
option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file"),
	make_option(c("-x", "--x_values"), type="character", 
		help="Name of column with values to be plotted"),
	make_option(c("-y", "--y_values"), type="character", 
		help="Name of column with values to be plotted"),
	make_option(c("-f", "--series"), type="character", 
		help="Name of column to be used as series"),
	make_option(c("-H", "--header"), action="store_false", default=TRUE,
        help="The input table not have header line"),
	 make_option(c("-X", "--x_title"), type="character", 
	 	help="Name of column to be used for bars titles"), 
	make_option(c("-Y", "--y_title"), type="character", 
	 	help="Title of y axis"),
	make_option(c("-F", "--output_format"), type="character", default="pdf", 
	 	help="pdf or jpeg file output format"),
        make_option(c("--x_max"), type="integer", default=NA,
                help="Use for set x axis max"),
        make_option(c("--x_min"), type="integer", default=NA,
                help="Use for set x axis min"),
        make_option(c("--y_max"), type="integer", default=NA,
                help="Use for set y axis max"),
        make_option(c("--y_min"), type="integer", default=NA,
                help="Use for set y axis min")
)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

data <- read.table(opt$data_file, sep="\t", header=opt$header)
if (opt$output_format == "pdf"){
	pdf(paste(opt$output, '.pdf', sep=""))
}else if(opt$output_format == "jpeg"){
	jpeg(paste(opt$output, '.jpeg', sep=""))
}	
	ggplot(data=data, aes(x=data[[opt$x_values]], y=data[[opt$y_values]], color=data[[opt$series]] ))  +
	geom_point() +
	xlim(opt$x_min, opt$x_max) +
	ylim(opt$y_min, opt$y_max) +
	xlab(opt$x_title) +
	ylab(opt$y_title) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_colour_discrete("Datasets")
dev.off()
