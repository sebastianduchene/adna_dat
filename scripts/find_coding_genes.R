library(ape)

coding_gene_info <- read.table('coding_genes_frames_moralis_ref_ancient.csv', row.names = 1, head = T, sep = ',', stringsAsFactors = F)


for(i in 1:nrow(coding_gene_info)){
      gene_temp <- coding_gene_info[i, 1]
      frame_temp <- coding_gene_info[i, 2]
      aln_temp <- read.dna(paste0('moralis_ref_ancient/', gene_temp), format = 'fasta')
      print(gene_temp)
      aln_frame <- aln_temp[, (frame_temp+1):ncol(aln_temp)]
      nearest_div_3 <- ncol(aln_frame) %% 3
      if(nearest_div_3 == 0){
        write.dna(aln_frame, file = paste0('coding_genes_moralis_ref_ancient/', gene_temp), nbcol = -1, colsep = '', format = 'fasta')
	write.dna(aln_frame, file = paste0('coding_genes_moralis_ref_ancient/', gsub('fasta', 'phy',gene_temp)), nbcol = -1, colsep = '')
      }else{
        aln_frame <- aln_frame[, 1:(ncol(aln_frame) - nearest_div_3)]
	write.dna(aln_frame, file = paste0('coding_genes_moralis_ref_ancient/', gene_temp), nbcol = -1, colsep = '', format = 'fasta')
	write.dna(aln_frame, file = paste0('coding_genes_moralis_ref_ancient/', gsub('fasta', 'phy',gene_temp)), nbcol = -1, colsep = '')
      }
}
