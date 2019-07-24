import argparse
import sys
import random
import operator
from itertools import imap


#5'-ACTGGCCGCTTCACTG-var-GGTACCTCTAGA-tag-AGATCGGAAGAGCGTCG-3'
#5'-F1-var-KpnI-filler_seq-XbaI-tag-R1-3'   total size = 230 nt

nucs = ['A', 'C', 'G', 'T']

KpnI = "GGTACC"
XbaI = "TCTAGA"

F1 = "ACTGGCCGCTTCACTG"
R1 = "AGATCGGAAGAGCGTCG"


# find hamming distasnce between two sequences
def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    ne = operator.ne
    return sum(imap(ne, seq1, seq2))


# filter once all oligos have been generated
def tag_generator(len_tags, num_tags):
    tags = []

    # Loop until we have 7500 tags that have hamming distance of at least 2
    while len(tags) < num_tags:
        tag = ''
        remove = False
        for i in range(len_tags):
            tag = tag + random.choice(nucs)

        for cur_tag in tags:
            if hamming_distance(tag, cur_tag) < 2:
                remove = True

        if not remove:
            tags.append(tag)

    return tags


def create_oligos(tags, str_permutations):
    pass

if __name__ == '__main__':
    all_tags = tag_generator(10, 7500)
    create_oligos(all_tags, )

