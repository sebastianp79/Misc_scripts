#!/bin/rscripts


#this script takes files with aggregate methylation data at loci an plots average line graphs as well as multiple individual graphs for example prepared with 'Aggregate_methylation_at_loci.R'
#input: aggregate data (columns (SampleID, bins, score)

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="location file name", metavar="character"),
  make_option(c("-p", "--out_path"), type="character", default=NULL,
              help="path the plots should be written to", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(dplyr)
library(ggplot2)


output_name <- paste(opt$out, "aggregated_Samples.pdf", sep = "_" )

#prepare aggregate plot

aggdata <- read.delim(file=opt$file, header=T, sep="\t")

aggdata %>%
        dplyr::group_by(bins) %>%
        dplyr::summarize(avgScore = mean(score, na.rm=T)) %>%
        ggplot(., aes(bins, avgScore)) +
           geom_line(size=1) +
           theme_bw() +
           xlab("distance in bp")  +
           ylab("average methylation level")

ggsave(output_name, path=opt$out_path, width=10, height=6)



#prepare plot with individual lines for each sample

output_name <- paste(opt$out, "individual_Samples.pdf", sep = "_" )

aggdata %>%
        ggplot(.,aes(bins, score, group=SampleID, color=SampleID)) +
           geom_line(size=1) +
           guides(color="none") +
           theme_bw() +
           xlab("distance in bp")  +
           ylab("average methylation level")
ggsave(output_name, path=opt$out_path, width=10, height=6)



#add heatmap

#get aggregate data and transform into matrix

library (gplots)

#number of samples found int that file
sample_num <- aggdata$SampleID %>%
				unique() %>%
				length()


mat <-  matrix(aggdata$score,
				nrow = sample_num,
				byrow=T)



names <- aggdata$SampleID %>%
				unique()

#plot the data using heatmap.2 and write to pdf in same directory as other output files
out_path_file <- paste(opt$out_path, opt$out, "_individual_samples_heatmap.pdf", sep = "" )

pdf(out_path_file)

heatmap.2(mat,
           Colv =FALSE,
           Rowv=TRUE,
          labRow=names,
           # data scaling
           scale = c("row"),
           na.rm=TRUE,
          trace=c("none"),
          col="bluered",
          labCol="none",
           density.info=c("none"))

dev.off()
