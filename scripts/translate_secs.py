import re, os
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


files = [i for i in os.listdir('moralis_ref_ancient') if 'fasta' in i]



coding_genes = list()
for f in files:
    print f
    dat_temp = Bio.SeqIO.to_dict(Bio.SeqIO.parse(open('moralis_ref_ancient/'+f, 'r'), 'fasta'))
    ref_temp = str(dat_temp['MethanobrevibacteroralisJMR01'].seq)
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
        continue


out_frame = pd.DataFrame(coding_genes)
out_frame.to_csv('coding_genes_frames_moralis_ref_ancient.csv')


