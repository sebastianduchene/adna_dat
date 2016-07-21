# Readme file

July 21 2016

- [This file](https://github.com/sebastianduchene/adna_dat/blob/master/summarise_coefVar.csv) contains the mean coefficient of rate variation for all genes with complete sequences for the six taxa. Those with a mean coefficients of rate variation < 1 were considered to be 'clocklike'.

- The genes with coefficient of rate variation < 1 were concatenated into [this alignment](https://raw.githubusercontent.com/sebastianduchene/adna_dat/master/concat_clock_genes.fasta)

- The concatenated alignment of clocklike genes was analysed in BEAST using a strict clock. To assess convergence we conducted two replicates of this anlalysis, with log files [here](https://github.com/sebastianduchene/adna_dat/blob/master/concat_clock_genes_replicate1.log) and [here](https://github.com/sebastianduchene/adna_dat/blob/master/concat_clock_genes_replicate2.log). The esimtates for the rate and age of the root node are the same. The rate has a mean of 1.92e-7 (HPD: 1.4341E-7, 2.4177E-7). The mean estimate for the age of the root is 7.09e5 (HPD: 5.4875E5, 9.215E5).

