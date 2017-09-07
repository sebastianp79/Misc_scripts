#!/bin/rscripts

#this script is meant to get aggregate methylation scores across loci of interest.

# it takes the locations of the loci of interests, a file containing the paths to methylation coverage files, the desired size of the bins and of the intervals

#use optparse for management of input arguments
library("optparse")

option_list = list(
  make_option(c("-l", "--file1"), type="character", default=NULL,
              help="location file name", metavar="character"),
  make_option(c("-m", "--file2"), type="character", default=NULL,
              help="path to methylation files", metavar="character"),
  make_option(c("-b", "--bin_size"), type="numeric", default=NULL,
              help="size of bins", metavar="numeric"),
  make_option(c("-r", "--region_size"), type="numeric", default=NULL,
              help="size of regions selected in file1", metavar="numeric"),
	make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#interrupt execution if file1 is not supplied (i shuld extend this to all files)
if (is.null(opt$file1)){
  print_help(opt_parser)
  stop("location file must be supplied.n", call.=FALSE)
}


#load the packages required to run this script
library("IRanges")
library("GenomicRanges")
library("dplyr")
library("genomation")
library("data.table")


#core function that takes in the methylation coverage GR, location GR, binsizes and regionsizes
# it returns a dafaframe
GenomicFeature.aggregateMeans <-  function(Methylation.GR, Location.GR, binsize=binsize, regionsize=regionsize){
  #Requires methylation file and location file as GenomicRanges objects

  bin.num=regionsize/binsize

  binMatrix <- Methylation.GR %>%
      ScoreMatrixBin(., Location.GR,
                      bin.num = bin.num,
                      bin.op = "mean",
                      strand.aware = TRUE,
                      weight.col = "propM",
                      is.noCovNA = TRUE)

  #prepare df with aggregate values for each bin position
  outdf <- data.frame(bins = seq(from=-(regionsize/2),
                                  to=(regionsize/2)-binsize,
                                  by=(regionsize/ncol(binMatrix))),
                        score = colMeans(binMatrix,na.rm=TRUE) )

  return(outdf)
}





#function to read in a file path of a methyation coverage file, import the file and return it as a GRanges object
meth.file.import <- function(meth_dir){

	#print out which file path is read right in </->now
	file_countr <- paste("Now reading file", basename(meth_dir), collapse=" ")
	print(file_countr)

  	Meth.file <- fread(input= meth_dir,
                  	data.table=FALSE,
                  	header=FALSE)

#convert file into GRanges object. use 1- starting coordinate
  Methylation.GR <- GRanges(seqnames = Meth.file$V1,
                          ranges = IRanges(start = Meth.file$V2,
                                          end = Meth.file$V3),
                          propM = Meth.file$V4)
rm(Meth.file)
return(Methylation.GR)
}










#wrapper function runs the functions meth.file.import and GenomicFeature.aggregateMeans and can be used with lapply (below)
get.sample.aggregates <- function(meth_dir,Location.GR, binsize= bin_size, regionsize=region_size ){

  #produce the GRanges object of methylation coverage file
  Methylation.GR <- meth.file.import(meth_dir)

  #run GenomicFeature.aggregateMeans to get coverage across a locus for a particular methylation file
  outdf <- GenomicFeature.aggregateMeans(Methylation.GR,
                                          Location.GR,
                                          binsize=5,
                                          regionsize=2000)
  return(outdf)

}




#import input files
#file containing the location
print(opt$file1)
Loc.file <- fread(input= opt$file1,
                  data.table=FALSE,
                  header=FALSE)

Location.GR <- GRanges(seqnames = Loc.file$V1,
                          ranges = IRanges(start = Loc.file$V2+1,
                                          end = Loc.file$V3),
                          strand = Loc.file$V6) %>%
                resize(., width=opt$r, fix="center")

#print the number of rows in the location file
print("number of rows in location file")
print(nrow(Loc.file))
rm(Loc.file)


#read file containing paths of files into R
path <- opt$file2
#this relies on the fact that the files are ending in .cov and that only the desired files are in the directory. That might not be ideal in the long run
paths <- dir(path, pattern = ".cov", full.names = TRUE)

#show which coverage files are used
print("Coverage files that will be used in this run")
paths

#run the commands and produce a single dataframe containing 3 columns (SampleID, bins, score)
outlist <- lapply(paths,
                  get.sample.aggregates,
                  Location.GR,
                  binsize=opt$b,
                  regionsize=opt$r)

names(outlist) <- basename(paths)

outlist.df <- rbindlist(outlist, idcol = "SampleID")

write.table(outlist.df,
            file=opt$out,
            row.names=F,
            col.names=T,
            quote=F,
            sep="\t")

warnings()

print(sessionInfo())
