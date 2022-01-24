#!/usr/bin/env Rscript

#
# make the dot plot from the coordinates produced by compare_kmers.py
#

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
inputFile = args[1]
outputFile = args[2]

df = read.table(inputFile,header=F,sep=" ",col.names=c("sequence1","sequence2"))

pdf(NULL)

ggplot(df, aes(x=sequence1,y=sequence2)) + geom_point()
ggsave(outputFile)
