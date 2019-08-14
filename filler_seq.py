import sys
import random


seq = ''
length = int(sys.argv[1])
nucs = ['A', 'C', 'G', 'T']
output_file = open(sys.argv[2], 'w')

# Generate a random sequence of length "length" and store in it an output file
for i in range(length):
    seq = seq + random.choice(nucs)

output_file.write(seq)
