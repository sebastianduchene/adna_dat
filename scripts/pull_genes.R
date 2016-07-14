library(ape)

#if(F){
# The alignment without gaps
aln <- read.dna('ungapped_aln_moralis.fasta', 'fasta')

dat <- read.table('indices_moralis.csv', head = T, sep = ',', row.names = 1)
annotations <- read.table('annotations_moralis.csv', head = T, row.names = 1, sep = ',')
#}

for(i in 3:nrow(annotations)){
      print(annotations[i, 1])
      start_i <- annotations[i, 3]
      end_i <- annotations[i, 4]
#      start_aln <- dat[dat[, 2] == start_i, 1]
#      end_aln <- dat[dat[, 2] == end_i, 1]
#      print(c(start_i, end_i, start_aln, end_aln))
      write.dna(aln[, start_i:end_i], file = paste0('ungapped_single_genes_moralis/',annotations[i, 1], '.fasta'), nbcol = -1, colsep = '', format = 'fasta')

}

