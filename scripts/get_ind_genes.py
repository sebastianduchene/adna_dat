import re, sys, os
import Bio
from Bio import SeqIO
import pandas as pd
import numpy as np


gff_data = open('GCF_000529525.1_Methanobrevibacter_oralis_JMR01_genomic.gff', 'r').readlines()

handle = open('Tree_MoralisMsmithii/Methanobrevibacter_progMauve_6_nocontigs.fasta', 'r')
complete_alignment = Bio.SeqIO.to_dict(Bio.SeqIO.parse(handle, 'fasta'))

ref_sequence = str(complete_alignment['MethanobrevibacteroralisJMR01'].seq)
print ref_sequence[:10]

annotation_data = list()
for l in gff_data:
    if '#' in l:
        continue
    if ('CDS' in l) and ('gene' in l):
        split_lines = re.split('\t', l)
        limits_ref = split_lines[3:5]
        gene_number = re.findall('gene[0-9]+', l)
        cds = re.findall('cds[0-9]', l)
        annotation_data.append([gene_number[0], cds[0], limits_ref[0], limits_ref[1]])

annotation_frame = pd.DataFrame(annotation_data)
annotation_frame.to_csv('annotations_moralis.csv')

index_aln = list()
counter = -1
for i in range(len(ref_sequence)):
    if ref_sequence[i] == '-':
        ref_index = 'null'
    else:
        ref_index = counter+1
        counter = ref_index

    index_aln.append([i, ref_index])

output_data = pd.DataFrame(index_aln)
output_data.to_csv('indices_moralis.csv')





