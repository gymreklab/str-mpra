import sys
import operator
from itertools import imap
from oligos import hamming_distance, reverse_complement


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
        assert len(motif) <= len(seq)
        for i in range(len(seq)-len(motif)+1):
            if hamming_distance(motif, seq[i:i+len(motif)]) < 2:
                print("Seq %s with subseq %s has too many similarities to motif %s"%(seq, seq[i:i+len(motif)], motif))
                break

