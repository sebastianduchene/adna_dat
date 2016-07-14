import re, os
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

files = [i for i in os.listdir('moralis_ref_ancient') if 'fasta' in i]

def write_dict(out_secs, out_file = 'test.fasta'):
    f = open(out_file, 'w')
    for i in out_secs:
        f.write('>'+i+'\n')
        f.write(out_secs[i]+'\n')
    f.close()
    
coding_genes = list()

for f in files:
#    print f
    dat_temp = Bio.SeqIO.to_dict(Bio.SeqIO.parse(open('moralis_ref_ancient/'+f, 'r'), 'fasta'))
    ref_temp = str(dat_temp['MethanobrevibacteroralisJMR01'].seq)
    dat_ancient = str(dat_temp['M_oralis_MAPCycle0_ConsensusStrictN_nocontigs'].seq)
    if ('N' in ref_temp) or ('n' in ref_temp):
        continue
    frame0 = str(Seq(ref_temp, IUPAC.unambiguous_dna).translate())[:-1]
    frame1 = str(Seq(ref_temp[1:], IUPAC.unambiguous_dna).translate())[:-1]
    frame2 = str(Seq(ref_temp[2:], IUPAC.unambiguous_dna).translate())[:-1]
    frame3 = str(Seq(ref_temp[3:], IUPAC.unambiguous_dna).translate())[:-1]
    frame4 = str(Seq(ref_temp[4:], IUPAC.unambiguous_dna).translate())[:-1]
    frame5 = str(Seq(ref_temp[5:], IUPAC.unambiguous_dna).translate())[:-1]
    check_stop_codons = [not('*' in i) for i in [frame0, frame1, frame2, frame3]]
    if any(check_stop_codons):
        correct_frame = np.where(check_stop_codons)[0][0]
        coding_genes.append([f, correct_frame])
        
    else:
        # Remove sites with stop codons for all frames and keep that with the most sites left
        find_stop_codons = [np.where([i == '*' for i in fr]) for fr in [frame0, frame1, frame2, frame3, frame4, frame5]]
        count_codons = np.array([len(i[0]) for i in find_stop_codons])
        min_codons = np.where(count_codons == min(count_codons))[0][0]
        min_codon_frame = [frame0, frame1, frame2, frame3, frame4, frame5][min_codons]
        
        rm_nucs = np.array([find_stop_codons[0][0]*3, find_stop_codons[0][0]*3+1,find_stop_codons[0][0]*3+2])
        rm_nucs = rm_nucs.flatten()
        
        nucs_no_stop_codons = ''.join([ref_temp[min_codons:][i] for i in range(len(ref_temp[min_codons:])) if not(i in rm_nucs)])
        nucs_no_stop_codons_ancient = ''.join([dat_ancient[min_codons:][i] for i in range(len(dat_ancient[min_codons:])) if not(i in rm_nucs)])
        
        out_secs = {'MethanobrevibacteroralisJMR01':nucs_no_stop_codons, 'M_oralis_MAPCycle0_ConsensusStrictN_nocontigs':nucs_no_stop_codons_ancient}
        write_dict(out_secs, 'coding_genes_moralis_ref_ancient/'+f)
        
