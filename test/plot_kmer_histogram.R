#!/usr/bin/env Rscript

#
# generally no plot will appear unless this is run manually from within R
#

library(ggplot2)

df = read.table("tmp/kmerz_test_data/10k_31mer_counts",header=F,sep=" ",col.names=c("kmer","rev","count"))

ggplot(df, aes(x=count)) +
  geom_histogram(binwidth=1)
  #scale_y_continuous(trans='log10')
