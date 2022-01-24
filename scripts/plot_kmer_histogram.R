#!/usr/bin/env Rscript

#
# plot kmer frequency histogram
#

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
inputFile = args[1]
outputFile = args[2]

df = read.table(inputFile,header=F,sep=" ",col.names=c("kmer","count"))

pdf(NULL)
ggplot(df, aes(x=count)) + geom_histogram(binwidth=1) +  scale_y_continuous(trans='log10')
ggsave(paste("log10_",outputFile,sep=""))

ggplot(df, aes(x=count)) + geom_histogram(binwidth=1)
ggsave(paste("linear_",outputFile,sep=""))
