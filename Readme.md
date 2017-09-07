

The script Aggregate_methylation_at_loci.R in combination with ggplot can be used to generate methylation density plots around genomic features for one or more files containing methylation data (or similar scores)


The script Aggregate_methylation_at_loci.R requires two input files: A bedfile containing locations of features (e.g. TSS, ChIP-seq peaks) in BED6 format (0-start), and a file containing the methylation data [chr, start, end, score] (1-start).

```{bash}
#Requires
#library("optparse")
#library("IRanges")
#library("GenomicRanges")
#library("dplyr")
#library("genomation")
#library("data.table")

Rscript --vanilla Aggregate_methylation_at_loci.R -l Feature.bed -m Dir/of/methylation_file -b bin_in_bp -r size_of_region -o Path/to/output/Aggregate_methylation.txt


```
The output file contains 3 columns: SampleID, bin (with distance to center in bp), methylation score. This output is formatted for easy plotting with ggplot


```{r}
library(ggplot)
library(dplyr)


#example for combine multiple samples
data <- read.delim(file="Path/to/output/Aggregate_methylation.txt", header = T, sep = "\t")
data %>%
        dplyr::group_by(bins) %>%
        dplyr::summarize(avgScore = mean(score, na.rm=T)) %>%
        ggplot(., aes(bins, avgScore)) +
           geom_line(size=1) +
           theme_bw() +
           xlab("distance in bp")  +
           ylab("average methylation level")
ggsave("Name-of-combined-aggregate-plot.pdf", path = "Path/to/plots", width = 7, height =5)



#example for individual sample (if only one sample was present in directory the above option will produce the same result)
data %>%
        filter(., SampleID == "SampleID") %>%
        ggplot(., aes(bins, score)) +
           geom_line(size=1) +
           theme_bw() +
           xlab("distance in bp")  +
           ylab("average methylation level")
ggsave("Name-of-individual-sample-aggregate-plotpdf", path = "Path/to/plots", width = 7, height =5)


```
