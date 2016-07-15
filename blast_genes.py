import os, re, sys
import Bio
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

q_seq = sys.argv[1]
seq_query = Bio.SeqIO.to_dict(Bio.SeqIO.parse(open(q_seq, 'r'), 'fasta'))
seq_ref = seq_query['MethanobrevibacteroralisJMR01']

result_handle = NCBIXML.read(NCBIWWW.qblast('blastn', 'nt', seq_ref.seq))

out_handle = open(re.sub('[.]fasta', 'blast_match.txt', q_seq), 'w')
out_handle.write(result_handle.alignments[0].title+'\n')
for i in result_handle.alignments[0].hsps:
    out_handle.write(i.query+'\n')
    out_handle.write(i.match+'\n')
    out_handle.write(i.sbjct+'\n')
out_handle.close()

