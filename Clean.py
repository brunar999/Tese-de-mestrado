Cleaning OTU fasta file 
from Bio import SeqIO
IN_FILE1 = ’97_OTUs_10.fasta’ IN_FILE2 = ’97_10_OTU_clean.txt’ OUT_FILE = ’out.fasta’

seq_ids = set()
with open(IN_FILE2) as fh:
for line in fh:
seq_id = line.split()[0] seq_ids.add(seq_id)
with open(OUT_FILE, ’w’) as oh:
for seq_record in SeqIO.parse(IN_FILE1, ’fasta’): if seq_record.id in seq_ids:
oh.write(seq_record.format(’fasta’))
