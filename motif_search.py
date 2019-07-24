import sys
import operator
from itertools import imap


seq_file = open(sys.argv[1], 'r')
comp_motif_file = open(sys.argv[2], 'r')
comp_motifs = []

seq = seq_file.readline().rstrip('\n')
for line in comp_motif_file:
    comp_motifs.append(line.rstrip('\n'))


def ham_dist(seq1, seq2):
    assert len(seq1) == len(seq2)
    ne = operator.ne
    return sum(imap(ne, seq1, seq2))


for motif in comp_motifs:
    if ham_dist(motif, seq) < 2:
        print("Seq has too many similarities to motif")
 
