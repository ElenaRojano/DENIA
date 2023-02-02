#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse, limma

################################################################
##                           CONFIGURE                        ##
################################################################
require(optparse)
require(limma)

option_list <- list(
  make_option(c("-n", "--input"), action="store", type="character",
              dest="input", help="Input file. Items files separated by commas. Files format: one column file without header"),
  make_option(c("-t", "--tag"), action="store", type="character", default = NULL,
              dest="tag", help="Sets tags sepparated by commas. If names are incorrect, auto names will be used"),
  make_option(c("-c", "--color"), action="store_true", type="logical", default = FALSE,
              dest="color", help="Flag to activate circle colours mode"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################
# Obtain sets file paths
sets_files <- unlist(strsplit(opt$input,","))
if(length(sets_files)<2){
	stop("There are less than 2 sets. Venn diagram needs 2 sets at least")
}

# Obtain sets tags
if(is.null(opt$tag)){
	sets_tags <- NULL
}else{
	sets_tags <- unlist(strsplit(opt$tag,","))
	if(length(sets_files) != length(sets_tags)){
		warning("Given tags have different dimensions than sets files given. Using regular tags")
		sets_tags <- NULL
	}	
}


# Load datasets
sets <- lapply(sets_files,function(file_path){
	message(paste("Loading set",file_path,"..."))
	set <- read.table(file = file_path, quote = "", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
	set <- as.vector(set[,1])
	return(set)
})






################################################################
##                          STUDY DATA                        ##
################################################################
# Obtain unique list of items
items <- unique(unlist(sets))

# Per each item, check sets
belonging <- unlist(lapply(sets, function(set){items %in% set}))
#write.table(belonging, file = "tabla")
#q()
belonging <- matrix(as.numeric(belonging),ncol = length(sets))

# Add set names
if(!is.null(sets_tags)){
	colnames(belonging) <- sets_tags
}
# rownames(belonging) <- items

# Calculate intersections
intersections <- vennCounts(belonging)
intersections[1,"Counts"] <- NA


################################################################
##                       PLOT AND EXPORT                      ##
################################################################
# Plot vennDiagram
pdf(paste(opt$output,"pdf",sep=".")) 
if(!opt$color){
	vennDiagram(intersections)
}else{
	vennDiagram(intersections,circle.col = rainbow(length(sets)))
}
dev.off()

