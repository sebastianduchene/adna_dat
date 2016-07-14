import re, os
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


files = [i for i in os.listdir('coding_genes') if 'fasta' in i]


aln_temp = Bio.SeqIO.to_dict(Bio.SeqIO.parse(open('coding_genes/'+files[0], 'r'), 'fasta'))

gaps_temp = list()
for s in aln_temp:
    s_translate = Seq(re.sub('-', 'A', str(aln_temp[s].seq))).translate()
    print s
    print s_translate



