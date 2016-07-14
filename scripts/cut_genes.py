import re, sys, os
import Bio
from Bio import SeqIO
import pandas as pd
import numpy as np

def ungap(aln, ref_sequence):
    #aln should be a dictionary
    non_gaps = np.where([i!='-' for i in ref_sequence])
    cut_sequence = dict()
    for s in aln:
        print 'converting '+s
        s_temp = str(aln[s].seq)
        s_cut = ''.join([s_temp[i] for i in non_gaps[0]])
        cut_sequence[s] = Bio.Seq.Seq(s_cut)
    return cut_sequence



gff_data = open('Methanobrevibacter_oralis/GCF_000529525.1_Methanobrevibacter_oralis_JMR01_genomic.gff', 'r').readlines()

handle = open('Methanobrevibacter_oralis/Tree_MoralisMsmithii/Methanobrevibacter_progMauve_6_nocontigs.fasta', 'r')
complete_alignment = Bio.SeqIO.to_dict(Bio.SeqIO.parse(handle, 'fasta'))

ref_sequence = str(complete_alignment['MethanobrevibacteroralisJMR01'].seq)

non_gaps = np.where([i=='-' for i in ref_sequence])

aln = complete_alignment

ungapped_aln = ungap(aln, ref_sequence)

out_handle = open('ungapped_aln_moralis.fasta', 'w')
for s in ungapped_aln:
    out_handle.write('>'+s+'\n')
    out_handle.write(str(ungapped_aln[s])+'\n')
out_handle.close()
