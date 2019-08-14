import sys
import re
from oligos import reverse_complement


seq_file = open(sys.argv[1], 'r')
initial_motifs = sys.argv[2].split(',')
comp_motifs = set()
seqs = []

for motif in initial_motifs:
    comp_motifs.add(reverse_complement(motif))
    comp_motifs.add(motif)

for line in seq_file:
    seqs.append(line.rstrip('\n'))

for seq in seqs:
    for motif in list(comp_motifs):
        if re.search(motif, seq):
            print("Sequence %s contains %s."%(seq, motif))

