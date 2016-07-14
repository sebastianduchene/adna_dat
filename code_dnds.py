import re, os
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def get_site_changes(seq_list):
    nuc_subs = 0
    nuc_total = 0
    for i in range(len(seq_list[0])):
        if not(('-' in [seq_list[0][i], seq_list[1][i]]) or ('n' in [seq_list[0][i], seq_list[1][i]])):
            if not(seq_list[0][i] == seq_list[1][i]):
                nuc_subs = nuc_subs+1
            nuc_total = nucs_total + 1
    if nuc_total == 0:
        return [0, 0, 0]
    else:
        return [nuc_subs, nuc_total, float(nuc_subs) / nuc_total]
    
def get_site_changes(seq_list):
    third = [i for i in range(2, len(seq_list[0]), 3)]
    first_second = [i for i in range(len(seq_list[0])) if i not in third]

    third_subs = 0
    for i in third:
        if any([not(q in [seq_list[0][i], seq_list[1][i]]) for q in ['-', 'n', 'x', 'X']]):
            if not(seq_list[0][i] == seq_list[1][i]):
                third_subs = third_subs+1
            
    first_second_subs = 0
    for i in first_second:
        if not(('-' in [seq_list[0][i], seq_list[1][i]]) or ('n' in [seq_list[0][i], seq_list[1][i]])):
            if not(seq_list[0][i] == seq_list[1][i]):
                first_second_subs = first_second_subs+1
    return [float(first_second_subs)*(2/3.), float(third_subs)*(1/3.)]
    
files = [i for i in os.listdir('coding_genes_moralis_ref_ancient/') if 'fasta' in i]

coding_genes_dnds = list()
for f in files:
    dat_temp = Bio.SeqIO.to_dict(Bio.SeqIO.parse(open('coding_genes_moralis_ref_ancient/'+f, 'r'), 'fasta'))
    seq_list = [str(dat_temp[i].seq) for i in dat_temp]
    nuc_div = get_site_changes(seq_list)
    if nuc_div[1] == 0:
        coding_genes_dnds.append([f, 0])
    else:
        coding_genes_dnds.append([f, nuc_div[0]/nuc_div[1]])
        
out_frame = pd.DataFrame(coding_genes_dnds)

out_frame.columns = ['gene_name', 'dn_ds']

out_frame.to_csv('coding_genes_moralis_ref_ancient_dnds.csv')