#!/usr/bin/env Rscript

#
# plot should appear as Rplot.pdf
#

library(ggplot2)

df = read.table("10k_31mer_counts",header=F,sep=" ",col.names=c("kmer","count"))

ggplot(df, aes(x=count)) +
  geom_histogram(binwidth=1)
  #scale_y_continuous(trans='log10')
