library(ape)

#Concatenate all genes with a mean coefficient of rate variation of less than 1
dat <- read.csv('summarise_coefVar.csv', head = T)
dat_clock <- dat[dat$meanCV < 1, ]

concat_data <- read.dna(gsub('_.+', '.fasta', dat$file_name[1]), format = 'fasta')
for(f in 2:nrow(dat_clock)){
      dat_temp <- read.dna(gsub('_.+', '.fasta', dat_clock$file_name[f]), format = 'fasta')
      concat_data <- cbind(concat_data, dat_temp)
      print(dat_temp)
}

write.dna(concat_data, file = 'concat_clock_genes.fasta', format = 'fasta', nbcol = -1, colsep ='')
